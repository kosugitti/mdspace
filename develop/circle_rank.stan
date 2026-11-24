// 円モデル(Okada & Imaizumi 1987)の順位尤度版(円モデル専用・試行用)
// ellipse_rank_std.stan の円モデル(use_size = 1・use_shape = 0)と同じ確率モデルで，母数の作りだけを変えた:
//   (1) 回転を止める: 計算の足場として点 1 を原点に置き，点 2〜N を cholesky_factor_cov[N-1, P](下三角・対角が正)で表す
//       (Stan Reference Manual "Cholesky factors of covariance matrices"。vonmises_mds.stan と同じ型)
//   (2) 円モデルで使わない母数(楕円の形 et・s_e)を置かない
//   (3) 半径は差だけを母数にする(2026-11-23): mu は全員の半径に同じ数を足しても変わらないので，
//       和が 0 の a_j を推定し，報告では OI と同じく r_j = a_j - min(a) と最小を 0 にずらす。
//       前の版(r = exp(c)・c = m_c + s_c * c_raw の対数正規)は半径の水準を事前分布だけで決めていて，
//       m_c と s_c の事後分布が尾根になった(相関 -0.92)
// 座標の事前分布は，前の版(正規化前の座標が独立に N(0, 1))から重心を積分で消したもの＝重心からのずれが N(0, 1)。
// 平行移動・回転・鏡映で変わらないので，足場の選び方は結果に効かない。
// 周辺尤度を前の版と同じ値にするため，回転と鏡映を積分で消した体積も残す(model ブロックの注)
// 潜在変数 zeta ~ N(mu, sig^2) を積分した順序プロビット(閾値 K-1 個・タイは一次の扱い)
// 閾値は cut = mean(mu) + sqrt(var(mu) + sig^2) * g，g は順序づけた N(0, 1)(閾値を潜在変数の広がりに合わせる)
// 足場の P+1 個の点は，互いに離れていて一つの超平面に乗らないものを選び，データの先頭に並べ替えて渡す
data{
    int<lower=2> N;                          // 対象の数
    int<lower=1, upper=N-1> P;               // 次元の数
    int<lower=1> L;                          // 観測セルの数(対角と欠測を除く)
    array[L] int<lower=1, upper=N> row;      // セルの行(評定する側 j)
    array[L] int<lower=1, upper=N> col;      // セルの列(評定される側 k)
    int<lower=2> K;                          // 非類似度の値の種類
    array[L] int<lower=1, upper=K> y;        // 各セルの非類似度が下から何番目の値か
}

parameters{
    cholesky_factor_cov[N-1, P] Lx;          // 点 2〜N の座標(正規化前・点 1 が原点)
    sum_to_zero_vector[N] a_raw;             // 半径の差(標準化・非中心化・和が 0)
    real<lower=0> s_a;                       // 半径の差の散らばり
    real<lower=0> sig;                       // 測定誤差(布置の広がりに比べた大きさ)
    ordered[K-1] g;                          // 閾値(潜在変数の広がりで標準化した単位)
}

transformed parameters{
    array[N, P] real tau;                    // 座標(重心が原点・平方和 N)
    vector[N] a = s_a * a_raw;               // 半径の差(和が 0)
    vector[L] mu;                            // モデルの値
    vector[K-1] cut;                         // 閾値(mu の単位)
    real ss = 0.0;                           // 正規化前の座標の，重心からのずれの平方和

    // 正規化: 重心を原点に，平方和を N に
    {
        array[P] real mn;
        for(p in 1:P){
            mn[p] = 0.0;
            for(j in 2:N){
                mn[p] = mn[p] + Lx[j-1, p];
            }
            mn[p] = mn[p] / N;
            ss = ss + mn[p]^2;               // 点 1(原点)のずれ
            for(j in 2:N){
                ss = ss + (Lx[j-1, p] - mn[p])^2;
            }
        }
        for(p in 1:P){
            tau[1, p] = -sqrt(N / ss) * mn[p];
            for(j in 2:N){
                tau[j, p] = sqrt(N / ss) * (Lx[j-1, p] - mn[p]);
            }
        }
    }

    // モデルの値 mu = delta - r_j + r_k(半径の水準は mu に効かないので，差 a で計算する)
    for(l in 1:L){
        real delta = 0.0;
        for(p in 1:P){
            delta = delta + (tau[row[l], p] - tau[col[l], p])^2;
        }
        delta = sqrt(delta);
        mu[l] = delta - a[row[l]] + a[col[l]];
    }

    // 閾値: 潜在変数 zeta = mu + 誤差 の平均と広がりに合わせる
    cut = mean(mu) + sqrt(variance(mu) + sig^2) * g;
}

model{
    // 座標の事前分布: 重心からのずれが N(0, 1)((N-1)P 次元の部分空間の上の正規分布)
    // 点 1 を原点に置いた座標からずれへの写像の体積が，軸ごとに N^(-1/2)
    target += -0.5 * (N - 1) * P * log(2 * pi()) - 0.5 * ss - 0.5 * P * log(N);
    // 回転と鏡映を積分で消した体積: (N-1) x P の行列 A = Lx Q(Q は直交行列)で dA = prod_p Lx[p,p]^(P-p) dLx dQ，
    // 直交群 O(P) の体積は 2^P pi^(P^2/2) / Gamma_P(P/2)(P = 2 なら 4 pi，P = 3 なら 16 pi^2)
    for(p in 1:P){
        target += (P - p) * log(Lx[p, p]);
    }
    target += P * log(2) + 0.5 * P^2 * log(pi()) - lmgamma(P, 0.5 * P);
    // 和が 0 のベクトルは N-1 次元の部分空間の上にあり，Stan の変換は等長(ヤコビアン 1)。N 成分に normal_lpdf を
    // かけると正規化定数が N 個ぶん付くので，余分な 1 個ぶん 0.5*log(2*pi) を戻す
    target += normal_lpdf(a_raw | 0, 1) + 0.5 * log(2 * pi());
    target += normal_lpdf(s_a | 0, 1) + log(2);   // 半正規
    target += normal_lpdf(sig | 0, 1) + log(2);
    target += normal_lpdf(g | 0, 1) + lgamma(K);  // 順序づけた K-1 個: (K-1)! を掛けて正規化

    // 尤度: zeta ~ N(mu, sig^2) を積分した順序プロビット
    target += ordered_probit_lpmf(y | mu / sig, cut / sig);
}

generated quantities{
    vector[N] r = a - min(a);                // 半径(OI と同じく最小を 0 にずらして報告する)
    vector[L] log_lik;                       // セルごとの対数尤度
    for(l in 1:L){
        log_lik[l] = ordered_probit_lpmf(y[l] | mu[l] / sig, cut / sig);
    }
}
