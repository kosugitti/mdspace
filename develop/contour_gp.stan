// 等高線マップモデル(Abelson 1954-55・小杉・藤原 2004)のガウス過程版(試作・2026-11-28)
// 布置は類似性データの MDS で決めたものを固定して渡す(集団 g ごとに別の布置でよい)。
// 人 i の対象 j への評定 v_ij を，布置の上の見えない好意の場の観測とみなす:
//   v_ij = alpha_i + f_g(x_gj) + h_i(x_gj) + e_ij,   e_ij ~ N(0, sig^2)
//   f_g ~ GP(0, eta0^2 k)   集団 g に共通の場
//   h_i ~ GP(0, eta_i^2 k)  人 i の場のずれ(h_i は積分で消し，人ごとに多変量正規)
//   k(d) = 1 / (1 + d^2 / (2 ell^2))   有理二次カーネルの alpha = 1(コーシー型)。
//   2 ell^2 = 1 で Abelson の 1 / (1 + d^2) と同じ
// 人ごとの整合の指標: 個人の場で説明できる割合 omega_i = eta_i^2 / (eta_i^2 + sig^2)
// 評定は当面は連続量として扱う(順序プロビットは次の段階)
data{
    int<lower=1> G;                          // 集団(布置)の数
    int<lower=2> J;                          // 対象の数
    array[G] matrix[J, 2] X;                 // 集団ごとの布置(2 次元)
    int<lower=1> N;                          // 人の数
    array[N] int<lower=1, upper=G> grp;      // 人の集団
    int<lower=1> L;                          // 観測した評定の数(欠測を除く)
    vector[L] v;                             // 評定(中央を 0 にしたもの)
    array[L] int<lower=1, upper=J> obj;      // 評定の対象
    array[N] int<lower=1> start;             // 人 i の評定が並びの何番目から始まるか
    array[N] int<lower=1> len;               // 人 i の評定の数
}

transformed data{
    array[G] matrix[J, J] D2;                // 対象の間の距離の二乗
    for(g in 1:G){
        for(j in 1:J){
            for(k in 1:J){
                D2[g, j, k] = squared_distance(X[g, j], X[g, k]);
            }
        }
    }
}

parameters{
    real<lower=0> ell;                       // 場の広がり(長さの尺度)
    real<lower=0> eta0;                      // 共通の場の大きさ
    array[G] vector[J] f_raw;                // 共通の場(非中心化)
    real m_eta;                              // 個人の場の大きさの対数の平均
    real<lower=0> s_eta;                     // 個人の場の大きさの対数の散らばり
    vector[N] eta_raw;                       // 個人の場の大きさ(非中心化)
    real<lower=0> sig;                       // 誤差
    real alpha0;                             // 評定の水準の平均
    real<lower=0> s_alpha;                   // 評定の水準の散らばり
    vector[N] alpha_raw;                     // 人ごとの評定の水準(非中心化)
}

transformed parameters{
    array[G] matrix[J, J] K;                 // カーネル行列(対角は 1)
    array[G] vector[J] f;                    // 共通の場(対象の点での値)
    vector[N] eta = exp(m_eta + s_eta * eta_raw);       // 個人の場の大きさ
    vector[N] alpha = alpha0 + s_alpha * alpha_raw;     // 人ごとの評定の水準
    for(g in 1:G){
        for(j in 1:J){
            for(k in 1:J){
                K[g, j, k] = 1.0 / (1.0 + D2[g, j, k] / (2.0 * square(ell)));
            }
        }
        f[g] = cholesky_decompose(add_diag(square(eta0) * K[g], 1e-8)) * f_raw[g];
    }
}

model{
    ell ~ lognormal(log(sqrt(0.5)), 0.5);   // 中央値は Abelson の値(2 ell^2 = 1)
    eta0 ~ normal(0, 1);
    m_eta ~ normal(-1, 1);
    s_eta ~ normal(0, 0.5);
    sig ~ normal(0, 1);
    alpha0 ~ normal(0, 1);
    s_alpha ~ normal(0, 1);
    for(g in 1:G){
        f_raw[g] ~ std_normal();
    }
    eta_raw ~ std_normal();
    alpha_raw ~ std_normal();
    for(i in 1:N){
        array[len[i]] int idx = obj[start[i]:(start[i] + len[i] - 1)];
        matrix[len[i], len[i]] S = add_diag(square(eta[i]) * K[grp[i]][idx, idx], square(sig));
        v[start[i]:(start[i] + len[i] - 1)] ~ multi_normal(alpha[i] + f[grp[i]][idx], S);
    }
}

generated quantities{
    vector[N] omega;                         // 個人の場で説明できる割合(整合の指標)
    vector[N] log_lik;                       // 人ごとの対数尤度(LOO 用)
    for(i in 1:N){
        array[len[i]] int idx = obj[start[i]:(start[i] + len[i] - 1)];
        matrix[len[i], len[i]] S = add_diag(square(eta[i]) * K[grp[i]][idx, idx], square(sig));
        omega[i] = square(eta[i]) / (square(eta[i]) + square(sig));
        log_lik[i] = multi_normal_lpdf(v[start[i]:(start[i] + len[i] - 1)] | alpha[i] + f[grp[i]][idx], S);
    }
}
