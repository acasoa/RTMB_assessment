
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
  logSigmaRec=0,
  logSigmaC=0,
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
  sigmaRec <- exp(logSigmaRec)
  sigmaC <- exp(logSigmaC)
  SL50 <- exp(logSL50)
  SL95 <- SL50 + exp(logdiff)
  
  ### length, weight and selex per age
  La <- linf*(1-exp(-k*(ages-t0)))
  Wa <- aw*La^bw
  Sa <- 1/(1+exp(-log(19)*((La-SL50)/(SL95-SL50))))
  
  ### Mortality matrixes
  Fmat <- outer(Sa,F)
  Zmat <- Fmat+M
  
  ### biomass matrix
  N <- matrix(0,nrow = Amax,ncol = length(years))
  
  # 1st row = recruitment
  N[1,1] <- R0
  for (a in 2:length(N[1,])) {
    N[1,a] <- R0*exp(RecDevs[a]-(0.5*sigmaRec^2))
  }
  
  # filling the matrix
  for (t in 1:(length(years)-1)) {
    for (a in 1:(Amax-1)) {
      N[a+1,t+1] <- N[a,t]*exp(-Zmat[a,t])
    }
    N[Amax,t+1] <- N[Amax-1,t]*exp(-Zmat[Amax-1,t]) + N[Amax, t] * exp(-Zmat[Amax, t])
  }
  
  ### predicted catch
  logPredC <- log(colSums(Fmat/Zmat*(1-exp(-Fmat-M))*N*Wa))
  
  -sum(dnorm(log(catch), logPredC, sigmaC, log=TRUE))
}

obj <- MakeADFun(nll,par, random = "RecDevs")
fit <- nlminb(obj$par,obj$fn,obj$gr)
fit$par

fit$convergence
max(abs(obj$gr(fit$par)))
sdrep$pdHess



