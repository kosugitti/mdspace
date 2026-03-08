library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
library(ltm)
library(smacof)
library(ggrepel)
# Purpose -----------------------------------------------------------------

## 1,双対尺度法の順序尺度
## 2.多次元展開法
## 3.GRMの確率モデル
## 4.GPCMの確率モデル
## 5.Ordinal MDS
## 6.Prefmap MDS
## 7.INDSCAL-FOLDING
## 8.状態空間MDS


# Dual Scaling ------------------------------------------------------------
Wisky <- matrix(c(1, 3, 6, 3, 5, 3, 6, 2, 0), nrow = 3)
DualScalingMA(Wisky)
DualScaling(Wisky)


# Dual Scaling for Ordinal Scale ------------------------------------------

ord.data <- matrix(c(
  3, 2, 3, 1, 2,
  6, 7, 5, 8, 7,
  8, 5, 7, 6, 4,
  1, 3, 2, 3, 5,
  4, 1, 1, 2, 1,
  5, 8, 6, 7, 8,
  2, 4, 4, 5, 3,
  7, 6, 8, 4, 6
), nrow = 5)
DualScalingOrdinal(ord.data)

source("http://aoki2.si.gunma-u.ac.jp/R/src/ro.dual.R", encoding = "euc-jp")
ro.dual(ord.data)


# Dual Scaling for Multiple choice ----------------------------------------


# Dual Scaling for paired comparison --------------------------------------


pair.data <- matrix(c(
  1, 0, 2, 1, 1, 2, 1, 1,
  2, 2, 2, 2, 2, 1, 2, 2,
  1, 2, 2, 0, 2, 1, 1, 2,
  2, 2, 2, 2, 2, 1, 1, 2,
  2, 1, 2, 1, 0, 2, 1, 2,
  1, 1, 1, 1, 1, 1, 1, 2
), ncol = 6)


# Multidimensional Unfolding model ---------------------------------------------------------

## simulation data
D <- 2
M <- 9
N <- 2000
Fsc <- matrix(nrow=N,ncol=D)
Apos <- matrix(nrow=M,ncol=D)
for(n in 1:N){
  for(d in 1:D){
    Fsc[n,d] = rnorm(1,0,1);
  }
}
Apos <- matrix(c(
                  0,-2,
                  0,-1,
                 -2, 0,
                  0, 0,
                  2, 0,
                  1, 1,
                 -2, 1,
                  0,-3,
                  1, 2),ncol=2,byrow=T)
plot(Apos[,1],Apos[,2])
beta0 <- 10
beta1 <- 1.5
sig <- 2

DataMat <- matrix(nrow=N,ncol=M)
for(n in 1:N){
  for(m in 1:M){
    tmp = 0
    for(d in 1:D){
      tmp = tmp + (Fsc[n,d]-Apos[m,d])^2
    }
    tmp = sqrt(tmp)
    DataMat[n,m] = beta0 - beta1*tmp + rnorm(1,mean=0,sd=sig)
  }
}

DataMat %>% as.data.frame %>% mutate(ID=row_number()) %>%
  tidyr::gather(key,val,-ID,factor_key=TRUE) %>%
  mutate(key=as.numeric(key)) -> dat.long

datalist <- list(
  L = NROW(dat.long), M = max(dat.long$key), N = max(dat.long$ID),
  Pid = dat.long$ID, Jid = dat.long$key, val = dat.long$val, D= 2
)

UNFLD <- stan_model("R/UNFLD.stan")
init.UNFLD <- vb(UNFLD,datalist)

init.UNFLD %>% rstan::extract() %>% data.frame %>%
  tidyr::gather(key,val) %>% group_by(key) %>%
  summarise(MAP=median(val)) %>% spread(key,MAP) %>%
  dplyr::select(starts_with('beta'),starts_with('a')) %>% as.data.frame -> init_list

init.mcmc <- list()
for(c in 1:4){
  init.mcmc[c] <- list(init_list)
}



### sign constraint
sign.const <- vb(UNFLD,datalist)
sign.const %>% rstan::extract() %>% data.frame %>%
  tidyr::gather(key,val) %>% group_by(key) %>%
  summarise(MED=median(val)) %>% spread(key,MED) %>%
  dplyr::select(starts_with('a')) %>% sign() %>% unlist %>% matrix(ncol=D,byrow=T) -> sign.list
s
datalist$SIGN <- sign.list
UNFLD.C <- stan_model("R/UNFLD_const.stan")
result.UNFLD <- sampling(UNFLD.C,datalist,iter=10000)
print(result.UNFLD,pars=c('a','beta0','beta1','sig'))
plot(result.UNFLD,pars=c('a','beta0','beta1','sig'),show_density=TRUE)
traceplot(result.UNFLD,pars=c('a','beta0','beta1','sig'))
result.UNFLD %>% rstan::extract() %>% data.frame() %>%
  dplyr::select(starts_with('a')) %>%
  tidyr::gather(key,val) %>% group_by(key) %>%
  summarise(MED=median(val),EAP=mean(val),MAP=bayestestR::map_estimate(val),
            L95 = quantile(val,0.025),U95=quantile(val,0.975)) %>%
  dplyr::select(MAP) %>% unlist() %>% matrix(ncol=D,byrow=T) %>%
  as.data.frame() %>% scale() %>% as.data.frame() %>%  rbind(Apos %>% as.data.frame() %>% scale %>% as.data.frame()) %>%
  dplyr::mutate(ID=c(LETTERS[1:9],letters[1:9])) %>%
  ggplot(aes(x=V1,y=V2,label=ID))+geom_point()+geom_text_repel()


# GRM ---------------------------------------------------------------------

GRMmodel <- stan_model("R/GRM.stan")



## Sample Data
dat <- Science[c(1, 3, 4, 7)]
result.grm <- grm(dat,IRT.param = TRUE)
dat %>%
  mutate(ID = row_number()) %>%
  tidyr::gather(key, val, -ID) %>%
  mutate(key = ifelse(key=="Comfort",1,
                      ifelse(key=="Work",2,
                             ifelse(key=="Future",3,4)))) %>%
  mutate(val = ifelse(val=="strongly disagree",1,
                      ifelse(val=="disagree",2,
                             ifelse(val=="agree",3,4))))  %>% arrange(ID)-> dat.long
datalist <- list(
  L = NROW(dat.long), C = max(dat.long$val), M = max(dat.long$key), N = max(dat.long$ID),
  Pid = dat.long$ID, Jid = dat.long$key, val = dat.long$val, D= 1
)

fit <- vb(GRMmodel,datalist)
result.grm
print(fit,pars=c("a","b","loc"))
fit %>% rstan::extract() %>% as.data.frame() %>%
  dplyr::select(starts_with("theta")) %>% tidyr::gather(key,val,factor_key=TRUE) %>%
  group_by(key) %>% summarise(EAP=mean(val)) %>%
  mutate(SC=factor.scores.grm(result.grm,resp.patterns = dat)$score.dat$z1) %>%
  ggplot(aes(x=EAP,y=SC))+geom_point()

## Sample Data2
datYG <- read_table("R/YG外向.txt") %>% rename(val=ABCDEFGHIJ) %>%
  mutate(A = str_sub(.$val,start=1, end=1) %>% as.numeric(),
         B = str_sub(.$val,start=2, end=2)%>% as.numeric(),
         C = str_sub(.$val,start=3, end=3)%>% as.numeric(),
         D = str_sub(.$val,start=4, end=4)%>% as.numeric(),
         E = str_sub(.$val,start=5, end=5)%>% as.numeric(),
         F = str_sub(.$val,start=6, end=6)%>% as.numeric(),
         G = str_sub(.$val,start=7, end=7)%>% as.numeric(),
         H = str_sub(.$val,start=8, end=8)%>% as.numeric(),
         I = str_sub(.$val,start=9, end=9)%>% as.numeric(),
         J = str_sub(.$val,start=10, end=10)%>% as.numeric()) %>%
  dplyr::select(-"val") %>%
  mutate(ID=row_number()) %>%
  tidyr::gather(key,val,-ID) %>%
  mutate(key=as.numeric(as.factor(key))) %>%
  mutate(val=val+1)

datalist <- list(
  L = NROW(datYG), C = max(datYG$val), M = max(datYG$key), N = max(datYG$ID),
  Pid = datYG$ID, Jid = datYG$key, val = datYG$val
)
fit <- sampling(GRMmodel,datalist)
print(fit,pars=c("a","loc"))

# GPCM --------------------------------------------------------------------
GPCM <- stan_model("R/GPCM.stan")
fit <- sampling(GPCM,datalist)
print(fit,pars=c("a","b","loc"))


# 多次元GRM ------------------------------------------------------------------
psych::bfi %>% dplyr::select(starts_with("C"),starts_with("E"),-education) %>%
  mutate(ID=row_number()) %>%
  # 逆転項目の処理
  mutate(C4 = 7-C4, C5 = 7-C5, E1 = 7-E1, E2=7-E2) %>%
  # 具合の良さそうな項目だけ残す
  dplyr::select(C1,C2,C3,E1,E2,ID) %>%
  # 数を減らす
  dplyr::slice(1:100) %>%
  tidyr::gather(key,val,-ID,factor_key=TRUE) %>% na.omit %>%
  mutate(key=as.numeric(key)) -> bfi.dat


#read_csv("R/CE.csv") %>% mutate(ID = row_number()) %>%
#  tidyr::gather(key,val,-ID,factor_key=TRUE) %>%
#  mutate(key = as.numeric(key),val=val+1) -> bfi.dat

datalist <- list(
  L = NROW(bfi.dat), C = max(bfi.dat$val), M = max(bfi.dat$key), N = max(bfi.dat$ID),
  Pid = bfi.dat$ID, Jid = bfi.dat$key, val = bfi.dat$val,D=2
)

init_list <- list(a=matrix(c(1.5,1.7,1.3,0.1,0.3,0.3,0.3,0.1,1.5,2.3),ncol=2))

GRMD <- stan_model("R/GRMD.stan")
fit_D <- sampling(GRMD,datalist,init=list(init_list,init_list,init_list,init_list))
print(fit_D,pars=c('a','b'))
# 綺麗なラベルスイッチング
traceplot(fit_D,pars=c("a[1,1]","a[1,2]"))
traceplot(fit_D,pars=c("a[2,1]","a[2,2]"))
traceplot(fit_D,pars=c("a[3,1]","a[3,2]"))
traceplot(fit_D,pars=c("a[4,1]","a[4,2]"))
traceplot(fit_D,pars=c("a[5,1]","a[5,2]"))

# 多次元GPCM -----------------------------------------------------------------
GPCMD <- stan_model("R/GPCMD.stan")
fit_pd <- sampling(GPCMD,datalist)
print(fit_pd,pars=c('a','b'))


# モデル比較 -------------------------------------------------------------------

## Sample Data
dat <- Science[c(1, 3, 4, 7)]
dat %>%
  mutate(ID = row_number()) %>%
  tidyr::gather(key, val, -ID) %>%
  mutate(key = ifelse(key=="Comfort",1,
                      ifelse(key=="Work",2,
                             ifelse(key=="Future",3,4)))) %>%
  mutate(val = ifelse(val=="strongly disagree",1,
                      ifelse(val=="disagree",2,
                             ifelse(val=="agree",3,4))))  %>% arrange(ID)-> dat.long
datalist <- list(
  L = NROW(dat.long), C = max(dat.long$val), M = max(dat.long$key), N = max(dat.long$ID),
  Pid = dat.long$ID, Jid = dat.long$key, val = dat.long$val, D= 2
)

dat.long %>% tidyr::spread(key,val) %>% dplyr::select(-ID) -> dat

GRM <- stan_model("R/GRM.stan")
GPCM <- stan_model("R/GPCM.stan")
UNFLD <- stan_model("R/UNFLD.stan")
source("R/DualScaling.R")
source("R/make_dominance.R")
source("R/DualScalingOrdinal.R")

result.DS <- DualScalingOrdinal(dat)
result.unfolding <- smacof::unfolding(dat,1)
result.grm <- grm(dat,IRT.param = TRUE)
result.GRM <- sampling(GRM,datalist)
result.GPCM <- sampling(GPCM,datalist)
result.UNFLD <- sampling(UNFLD,datalist,init=init.mcmc)

## GRMは多次元にしたらラベルスイッチング
## UNFLDは制約のかけ方に問題？収束が悪い。vbだと初期値に依存してる？答えが毎回違う。

result.DS$NormedCol
result.grm
result.unfolding
result.unfolding$conf.col
result.GRM %>% print(pars=c('a','b'))
result.GPCM %>% print(pars=c('a','b'))
result.UNFLD %>% print(pars=c('a','beta0','beta1','phi[1,3]'))
