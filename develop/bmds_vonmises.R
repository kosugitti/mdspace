rm(list=ls())
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(circular)
# circluarパッケージの使い方になれる
vec <- c()
for(d in seq(-pi,pi,0.1)){
  vec <- c(vec,dvonmises(circular(d),mu=circular(0),kappa=circular(1)))
}
plot(vec)


# 推定モデルを確認する --------------------------------------------------------------

# 各点の平均と分散
psi <- c(-3,1.5,0.8,3)
phi <- c(0,1,2,3)

bins <- 200
Angs <- seq(-pi,pi,length.out = bins)

df <- data.frame(ID=rep(1:4,each=bins))

dens <- c()
for(i in 1:4){
  for(k in 1:length(Angs))
  ## 角度を求める
  ## tan^{-1} = \frac{y2-y1}{x2-x1}
  dens <- c(dens,dvonmises(circular(Angs[k]),circular(psi[i]),circular(phi[i])))
}

sig <- 0.3

df$ang <- rep(Angs,4)
df$dens <- dens
df$dat <- dens + rnorm(4*bins,0,sig)

model <- stan_model("develop/vMtest2.stan")
standata <- list(N=NROW(df),P=4,Idx=df$ID,Ang=df$ang,D=df$dat)
fit <- sampling(model,standata)
print(fit,pars=c("mu","kappa","sig"))
