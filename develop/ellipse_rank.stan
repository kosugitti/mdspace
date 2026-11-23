// 楕円モデル(Okada 1990 の Model IIB)の順位尤度版
// 潜在変数 zeta ~ N(mu, sig^2) を積分した順序プロビット(閾値 K-1 個・タイは一次の扱い)
// mu_jk = delta_jk - rho_jk + rho_kj の符号は縛らない(負＝j の楕円が j-k を結ぶ直線の上で k を呑み込む配置も許す)
// (2026-11-21: mu > 0 を母数の空間で課すと，職業データでは事後分布の約 98% が制約の外にあり全遷移が発散した。
//  制約を外すと発散 1%・Rhat 1.02 以下。P(mu_jk < 0) をセルごとに報告して非対称が対称部を超える組を診断する)
// 布置は平方和 N に正規化して sig を推定し，軸ごとに基準の対象の座標を正に組み込む(鏡映を止める)
// (2026-11-21: 対数の版 log(zeta) ~ N(log(mu), 1) は最小のセルの mu が 0 へ吸い寄せられて全遷移が発散した)
// 楕円は log(phi_jp) = c_j + s_e * et_jp(大きさ c_j と，面積を保つ形のずれ et_j)
// use_size・use_shape で下位モデルを選ぶ: (1,1) 楕円モデル，(1,0) 円モデル，(0,1) 等面積の楕円モデル，(0,0) 対称
// (使わない側の散らばり s_c・s_e は事前分布だけで動く)
data{
    int<lower=2> N;                          // 対象の数
    int<lower=1> P;                          // 次元の数
    int<lower=1> L;                          // 観測セルの数(対角と欠測を除く)
    array[L] int<lower=1, upper=N> row;      // セルの行(評定する側 j)
    array[L] int<lower=1, upper=N> col;      // セルの列(評定される側 k)
    int<lower=2> K;                          // 非類似度の値の種類
    array[L] int<lower=1, upper=K> y;        // 各セルの非類似度が下から何番目の値か
    array[P] int<lower=1, upper=N> anchor;   // 軸ごとの基準の対象(その軸で重心から遠いもの)
    int<lower=0, upper=1> use_size;          // 1: 大きさが対象ごとに違う / 0: 全員同じ大きさ
    int<lower=0, upper=1> use_shape;         // 1: 楕円 / 0: 円
}

parameters{
    array[N-1, P] real yr;                   // 基準以外の対象の座標(正規化前)
    array[P] real v;                         // 基準の対象の位置(正規化前・他の平均 + exp(v))
    array[N] real c_raw;                     // 楕円の大きさ(標準化・非中心化)
    array[N] sum_to_zero_vector[P] et;       // 楕円の形(標準化したずれ・和が 0)
    real m_c;                                // 大きさの平均
    real<lower=0> s_c;                       // 大きさの散らばり
    real<lower=0> s_e;                       // 形の散らばり(0 なら円モデル)
    real<lower=0> sig;                       // 測定誤差(布置の広がりに比べた大きさ)
    ordered[K-1] gam;                        // 閾値(sig の単位)
}

transformed parameters{
    array[N, P] real tau_s;                  // 正規化前の座標
    array[N, P] real tau;                    // 座標(重心が原点・平方和 N)
    array[N, P] real<lower=0> phi;           // 半軸
    array[N] real c;                         // 楕円の大きさ(半軸の幾何平均の対数)
    vector[L] mu;                            // モデルの値

    // 正規化前の座標: 基準の対象は，他の対象の平均 + exp(v)
    for(p in 1:P){
        real tmp = 0.0;
        int i = 0;
        for(j in 1:N){
            if(j != anchor[p]){
                i = i + 1;
                tau_s[j, p] = yr[i, p];
                tmp = tmp + yr[i, p];
            }
        }
        tau_s[anchor[p], p] = tmp / (N - 1) + exp(v[p]);
    }

    // 正規化: 重心を原点に，平方和を N に
    {
        real ss = 0.0;
        array[P] real mn;
        for(p in 1:P){
            mn[p] = 0.0;
            for(j in 1:N){
                mn[p] = mn[p] + tau_s[j, p];
            }
            mn[p] = mn[p] / N;
            for(j in 1:N){
                ss = ss + (tau_s[j, p] - mn[p])^2;
            }
        }
        for(p in 1:P){
            for(j in 1:N){
                tau[j, p] = sqrt(N / ss) * (tau_s[j, p] - mn[p]);
            }
        }
    }

    // 半軸: 大きさ c_j と形 s_e * et_jp(どちらも非中心化。中心化の c ~ normal(m_c, s_c) は
    // 事後確率最大値で s_c -> 0 へ潰れる(事前密度が無限大になる)ので使わない)
    for(j in 1:N){
        c[j] = m_c + use_size * s_c * c_raw[j];
    }
    for(j in 1:N){
        for(p in 1:P){
            phi[j, p] = exp(c[j] + use_shape * s_e * et[j][p]);
        }
    }

    // モデルの値 mu = delta - rho_jk + rho_kj
    for(l in 1:L){
        real delta = 0.0;
        real tj = 0.0;
        real tk = 0.0;
        for(p in 1:P){
            real dif = tau[row[l], p] - tau[col[l], p];
            delta = delta + dif^2;
            tj = tj + (dif / phi[row[l], p])^2;  // 楕円 j の k の方向
            tk = tk + (dif / phi[col[l], p])^2;  // 楕円 k の j の方向
        }
        delta = sqrt(delta);
        mu[l] = delta - delta / sqrt(tj) + delta / sqrt(tk);
    }
}

model{
    // 事前分布(周辺尤度を計算するので，~ ではなく target += *_lpdf で正規化定数を残す。
    // 2026-11-22: ~ は -0.5*log(2*pi) や log(gam_sd) を落とし，次元を 1 つ増やすと約 18 が引かれずに残って
    // 周辺尤度が次元とともに増え続けた)
    target += normal_lpdf(to_array_1d(yr) | 0, 1);
    for(p in 1:P){
        target += normal_lpdf(tau_s[anchor[p], p] | 0, 1);  // 基準の対象も同じ事前分布
    }
    target += sum(v);                        // tau_s[anchor] = 平均 + exp(v) のヤコビアン
    target += normal_lpdf(c_raw | 0, 1);
    for(j in 1:N){
        // 和が 0 のベクトルは P-1 次元の部分空間の上にあり，Stan の変換は等長(ヤコビアン 1)。P 成分に normal_lpdf を
        // かけると正規化定数が P 個ぶん付くので，余分な 1 個ぶん 0.5*log(2*pi) を戻す(2026-11-22)
        target += normal_lpdf(et[j] | 0, 1) + 0.5 * log(2 * pi());
    }
    target += normal_lpdf(m_c | 0, 1);
    target += normal_lpdf(s_c | 0, 1) + log(2);   // 半正規
    target += normal_lpdf(s_e | 0, 1) + log(2);
    target += normal_lpdf(sig | 0, 1) + log(2);
    target += normal_lpdf(gam | 0, 5) + lgamma(K);  // 順序づけた K-1 個: (K-1)! を掛けて正規化

    // 尤度: zeta ~ N(mu, sig^2) を積分した順序プロビット
    target += ordered_probit_lpmf(y | mu / sig, gam);
}

generated quantities{
    vector[L] log_lik;                       // セルごとの対数尤度(LOO 用)
    for(l in 1:L){
        log_lik[l] = ordered_probit_lpmf(y[l] | mu[l] / sig, gam);
    }
}
