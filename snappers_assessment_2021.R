library(RTMB)

dat <- read.table("ex0.dat", header = TRUE)

par <- list(
  logN1=log(500),
  logF=log(0.2),
  logsigma=0
)

nll <- function(par){
  getAll(par,dat)
  N1 <- exp(logN1)
  F <- exp(logF)
  M <- 0.2
  sigma <- exp(logsigma)
  
  N <- c(N1 , rep (0 ,length(catch)-1))
  for (i in 2:length(N)) {
    N[i] <- N[i-1]*exp(-F-M)
  }
  
  logPredC <- c(rep (0 ,length(catch)))
  for (i in 1:length(logPredC)) {
    logPredC[i] <- log(F/(F+M)*(1-exp(-F-M))*N[i])
  }
  
  ADREPORT(N)
  
  -sum(dnorm(x=log(catch),mean=logPredC,sd=sigma,log = TRUE))
}

obj <- MakeADFun(nll, par)
fit <- nlminb(obj$par,obj$fn,obj$gr)
exp(fit$par)
sdrep <- sdreport(obj)
pl <- as.list(sdrep, "Est")
plr <- as.list(sdrep, "Est", report=TRUE)

######
#data
cat <- read.csv("data/analis.sc1.csv")
ind <- read.csv("data/analis.index.csv")
len <- read.csv("data/analis.LenComp.csv", header = F)
colnames(len) <- as.character(len[1,])
len <- len[-1,]
inputs <- readRDS(paste0(getwd(),"/data/input.rds"))
cat
ind
len
inputs


dat <- list(
  years=as.vector(cat$Year),
  catch=as.vector(cat$Handline),
  index.obs=as.vector(ind$Index),
  index.yrs=as.vector(ind$Year),
  index.idx=as.vector(match(ind$Year, cat$Year)),
  length.obs=as.matrix((as.matrix(len[, 6:ncol(len)]))),
  length.yrs=as.vector(as.numeric(len$Year)),
  length.idx=as.vector(match(len$Year,cat$Year)),
  length.nsamps=as.vector(as.numeric(len$Nsamps)),
  length.bins=as.vector(as.numeric(colnames(len)[6:ncol(len)])),
  linf=inputs$Linf_f_fix,
  k=inputs$k_f_fix,
  t0=inputs$t0_f_fix,
  aw=inputs$WLa_f_fix,
  bw=inputs$WLb_f_fix,
  M=inputs$M_f_fix,
  Amax=30,
  ages=as.vector(1:30),
  SL50=inputs$L50_f_fix,SL95=inputs$L95_f_fix,
  CV_len = as.vector(inputs$CV_lt_f_old_mean),
  sigma_R = inputs$sigmaR
)

par <- list(
  logR0=log(10000),
  logq=log(0.001),
  logF=log(0.2),
  logSigmaIdx=0,
  logSigmaRec=0,
  RecDev=rep(0,length(dat$years))
)


#### model


nll <- function(par){
  getAll(par,dat)
  R0 <- exp(logR0)
  q <- exp(logq)
  F <- exp(logF)
  sigmaIdx <- exp(logSigmaIdx)
  sigmaRec <- exp(logSigmaRec)
  
  N <- matrix(0,nrow = Amax,ncol = length(years))
  N[1,1] <- R0
  for (i in 2:length(N[1,])) {
    N[1,i] <- R0*exp(RecDev[i]-(0.5*sigmaRec^2))
  }
  
}

































