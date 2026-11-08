
#####data
cat <- read.csv("data/analis.sc1.csv")
ind <- read.csv("data/analis.index.csv")
len <- read.csv("data/analis.LenComp.csv", header = F)
colnames(len) <- as.character(len[1,])
len <- len[-1,]
inputs <- readRDS(paste0(getwd(),"/data/input.rds"))

##### input data
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
  CV_len = as.vector(inputs$CV_lt_f_old_mean),
  sigma_R = inputs$sigmaR
)

##### parameters list
par <- list(
  logR0=log(100000),
  logq=log(0.001),
  logF=rep(log(0.2),length(dat$years)),
  logSigmaIdx=0,
  #logSigmaRec=0,
  #logSigmaC=0,
  RecDevs=rep(0,length(dat$years)),
  logSL50=log(inputs$L50_f_fix),
  logdiff=log(inputs$L95_f_fix-inputs$L50_f_fix)
)

##### model
nll <- function(par){
  getAll(par,dat)
  R0 <- exp(logR0)
  q <- exp(logq)
  F <- exp(logF)
  sigmaIdx <- exp(logSigmaIdx)
  sigmaRec <- 0.6
  sigmaC <- 0.1
  SL50 <- exp(logSL50)
  SL95 <- SL50 + exp(logdiff)
  
  ### length, weight and selex per age
  La <- linf*(1-exp(-k*(ages-t0)))
  Wa <- aw*La^bw/1e6
  Sa <- 1/(1+exp(-log(19)*((La-SL50)/(SL95-SL50))))
  
  ### Mortality matrixes
  Fmat <- outer(Sa,F)
  Zmat <- Fmat+M
  
  ### biomass matrix
  N <- matrix(0,nrow = Amax,ncol = length(years))
  
  ret <- 0
  
  # 1st row = recruitment
  N[1,1] <- R0
  for (t in 2:length(N[1,])) {
    N[1,t] <- R0*exp(RecDevs[t]-(0.5*sigmaRec^2))
  }
  
  ret <- ret - sum(dnorm(RecDevs,0,sigmaRec,log=TRUE))
  
  # filling the matrix
  for (t in 1:(length(years)-1)) {
    for (a in 1:(Amax-1)) {
      N[a+1,t+1] <- N[a,t]*exp(-Zmat[a,t])
    }
    N[Amax,t+1] <- N[Amax-1,t]*exp(-Zmat[Amax-1,t]) + N[Amax, t] * exp(-Zmat[Amax, t])
  }
  
  ### predicted catch
  logPredC <- log(colSums(Fmat/Zmat*(1-exp(-Fmat-M))*N*Wa))
  
  ret <- ret - sum(dnorm(log(catch), logPredC, sigmaC, log=TRUE))
  
  ### abundance index
  Ba <- N*Wa*Sa
  IndexPred <- colSums(q*Ba)
  IndexPred <- IndexPred[index.idx]
  
  ret <- ret - sum(dnorm(log(index.obs), log(IndexPred), sigmaIdx, log=TRUE))
  
  ### lengths
  low <- c(-Inf,rep(0,length(length.bins)-1))
  for (i in 2:length(length.bins)) {
    low[i] <- length.bins[i]-0.5
  }
  
  upp <- c(rep(0,length(length.bins)-1), Inf)
  for (i in 1:(length(length.bins)-1)) {
    upp[i] <- length.bins[i]+0.5
  }
  
  Lmat <- matrix(0,nrow=length(ages), ncol=length(length.bins))
  for (t in 1:length(length.bins)) {
    for (a in 1:Amax) {
      Lmat[a,t] <-  pnorm((Lmat[a,t]+upp[t]),La[a],CV_len*La[a]) - pnorm((Lmat[a,t]+low[t]),La[a],CV_len*La[a])
    }
  }
  
  Lmat <- Lmat/rowSums(Lmat)
  
  Cnum <- (Fmat/Zmat)*(1-exp(-Zmat))*N
  Pmat <- matrix(0,nrow=length(years), ncol=length(length.bins))
  for (t in 1:length(years)) {
    Lvul <- colSums(Cnum[,t]*Lmat)
    Lvul <- Lvul/sum(Lvul)
    Pmat[t,] <- Lvul
  }
  
  Pmat <- Pmat[length.idx,]
  
  ret_len <- rep(0,nrow(Pmat))
  for (i in 1:length(ret_len)) {
    ret_len[i] <- -sum(dmultinom(x=as.vector(length.obs[i,]), prob=Pmat[i,], log=TRUE))
  }
  
  ret <- ret + sum(ret_len)
  
  ### reports
  B=colSums(N*Wa)
  ADREPORT(B)
  ADREPORT(F)
  ADREPORT(R0)
  
  return(ret)
}

obj <- MakeADFun(nll,par, random = "RecDevs")
fit <- nlminb(obj$par,obj$fn,obj$gr)
fit$par

fit$convergence
max(abs(obj$gr(fit$par)))
sdrep$pdHess

sdrep <- sdreport(obj)
pl <- as.list(sdrep, "Est")
plr <- as.list(sdrep, "Est", report=TRUE)


