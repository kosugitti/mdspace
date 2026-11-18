// フォンミーゼスモデル(AMISESCAL; 荘島 2011, 2012)のベイズ版 ----------------------------
// 設計メモ: develop/設計メモ_フォンミーゼスモデル_20261118.md
//
// モデル(2次元):
//   theta_ij = x_i から x_j を見た角度, delta_ij = ||x_i - x_j||
//   pi_ij    = c_i * exp{kappa_i (cos(theta_ij - mu_i) - 1)}   (c_i = 山の高さ = pi の最大値)
//   xi_ij    = (1 - pi_ij) * delta_ij
//   g_ij     ~ normal(xi_ij, sigma)
// omega_fixed = 1 なら荘島 (2011) の通常の vM 密度(omega = 1/(2 pi)):
//   c_i = exp(kappa_i) / (2 pi I0(kappa_i)) で kappa_i に従属し，c_i < 1 のため kappa_i < 6.5508
// omega_fixed = 0 なら c_i は自由(荘島 2012 の omega_i)。c_i ~ beta で退化を抑える
//
// 識別(回転・鏡映・平行移動): 計算の足場として点1を原点に置き，残りの N-1 点を cholesky_factor_cov[N-1, 2]
// (下三角・対角が正)で表す。点2は第1軸の正の側，点3は第2軸の正の側に乗る
// (Stan Reference Manual "Cholesky factors of covariance matrices")。
// mu は布置と一緒に回るので，布置を縛れば mu も識別される
// 基準の3点は互いに離れていて一直線に並ばないものを選び，データの先頭に並べ替えて渡す

data {
  int<lower=3> N;                       // 対象の数
  int<lower=1> L;                       // 観測したセルの数(対角と欠測を除く)
  array[L] int<lower=1, upper=N> row;   // 行(i: 評価する側)
  array[L] int<lower=1, upper=N> col;   // 列(j: 評価される側)
  vector[L] g;                          // 非類似度(尺度をそろえたもの)
  int<lower=0, upper=1> omega_fixed;
  real<lower=0> kappa_max;              // omega_fixed = 1 なら 6.5508 未満にする
  real<lower=0> prior_x_sd;             // 座標の事前分布の SD
  real<lower=0> prior_c_a;              // c ~ beta(a, b)(omega_fixed = 0 のとき)
  real<lower=0> prior_c_b;
}

parameters {
  cholesky_factor_cov[N - 1, 2] Lx;     // 点2〜N の座標(下三角・対角が正)
  array[N] unit_vector[2] m;            // 平均方向(単位ベクトル)
  vector<lower=0, upper=kappa_max>[N] kappa;
  vector<lower=0, upper=1>[omega_fixed ? 0 : N] c_free;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[N, 2] X;
  vector[N] mu;
  vector[N] c;
  vector[L] xi;
  X[1] = rep_row_vector(0, 2);
  X[2:N] = Lx;
  for (n in 1:N) {
    mu[n] = atan2(m[n][2], m[n][1]);
  }
  if (omega_fixed) {
    // exp(kappa) / (2 pi I0(kappa)) を桁あふれしないよう log で
    for (n in 1:N) {
      c[n] = exp(kappa[n] - log(2 * pi()) - log_modified_bessel_first_kind(0, kappa[n]));
    }
  } else {
    c = c_free;
  }
  for (l in 1:L) {
    row_vector[2] d = X[col[l]] - X[row[l]];
    real delta = sqrt(dot_self(d) + 1e-12);
    // cos(theta - mu) = (d / delta) . m
    real cs = (d[1] * m[row[l]][1] + d[2] * m[row[l]][2]) / delta;
    real p = c[row[l]] * exp(kappa[row[l]] * (cs - 1));
    xi[l] = (1 - p) * delta;
  }
}

model {
  // 事前分布
  // 事前分布は重心からのずれに置く: 平行移動・回転・鏡映で変わらないので，
  // 点1を原点に置く縛りは計算の足場にすぎず，どの点を足場にしても結果は同じ
  {
    row_vector[2] xbar = rep_row_vector(1.0 / N, N) * X;
    to_vector(X - rep_matrix(xbar, N)) ~ normal(0, prior_x_sd);
  }
  // 回転の縛りの補正: 点2を第1軸に乗せるのは回転角 phi を積分で消すことにあたり，
  // 極座標の体積要素 a da dphi から a(点1から点2までの距離)が一つ出る。
  // これを足すと事後分布が足場の選び方によらなくなる(平行移動と鏡映には出ない)
  target += log(Lx[1, 1]);
  kappa ~ gamma(2, 0.5);                // 平均4・裾は kappa_max で切る
  if (!omega_fixed) {
    c_free ~ beta(prior_c_a, prior_c_b);
  }
  sigma ~ student_t(4, 0, 1);
  // 尤度
  g ~ normal(xi, sigma);
}

generated quantities {
  // 推定の後で重心を原点に移した布置(足場の点1を特別扱いしない)
  matrix[N, 2] Xc = X - rep_matrix(rep_row_vector(1.0 / N, N) * X, N);
  vector[L] log_lik;
  for (l in 1:L) {
    log_lik[l] = normal_lpdf(g[l] | xi[l], sigma);
  }
}
