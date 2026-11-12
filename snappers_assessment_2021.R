setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(RTMB)
library(SPoRC)
library(scales)
library(dplyr)
#### data preparation ####
cat <- read.csv("data/analis.sc1.csv")
ind <- read.csv("data/analis.index.csv")
len <- read.csv("data/analis.LenComp.csv", header = FALSE)
colnames(len) <- as.character(len[1,])
len <- len[-1,]
inputs <- readRDS(paste0(getwd(),"/data/input.rds"))

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
  bw=inputs$WLb_f_fix
  #,M=inputs$M_f_fix
  ,Amax=inputs$Nages_in,
  ages=as.vector(1:inputs$Nages_in),
  CV_len = as.vector(inputs$CV_lt_f_old_mean),
  sigma_R = inputs$sigmaR,
  L50=inputs$L50_f_fix,
  L95=inputs$L95_f_fix
  #,h=inputs$h
  ,rec=1,
  ESS=FALSE
)

par <- list(
  logR0=log(10000000),
  logq=log(0.001),
  logF=rep(log(0.2),length(dat$years)),
  logSigmaIdx=0,
  #logSigmaRec=0,
  #logSigmaC=0,
  RecDevs=rep(0,length(dat$years)),
  logSL50=log(as.numeric(inputs$Sel50)),
  logdiff=log(as.numeric(inputs$Selpeak)-as.numeric(inputs$Sel50)),
  h=inputs$h,
  M=inputs$M_f_fix
)

#### model configuration ####
nll <- function(par){
  getAll(par,dat)
  ### loading data ####
  if (!exists("rec")) {
    warning("Recruitment assumed flat", call. = FALSE)
    rec <- 0
  }else if (rec==1) {
    warning("Recruitment assumed Beverton-Holt", call. = FALSE)
  } else if (rec==2) {
    warning("Recruitment assumed Ricker", call. = FALSE)
  }
  R0 <- exp(logR0)
  q <- exp(logq)
  Fi <- exp(logF)
  sigmaIdx <- exp(logSigmaIdx)
  sigmaRec <- sigma_R
  sigmaC <- 0.2
  SL50 <- exp(logSL50)
  SL95 <- SL50 + exp(logdiff)
  
  
  ### length, weigth and selex per age ####
  La <- linf*(1-exp(-k*(ages-t0)))
  Wa <- aw*La^bw/1e6
  Sa <- 1/(1+exp(-log(19)*((La-SL50)/(SL95-SL50))))
  Mata = 1.0 / (1.0 + exp(-log(19.0) * (La - L50) / (L95 - L50)))
  
  ### Mortality matrixes ####
  Fmat <- outer(Sa,Fi)
  Zmat <- Fmat+M
  
  ### biomass matrix ####
  N <- matrix(0,nrow = Amax,ncol = length(years))
  
  ### ssb virgin ####
  SSBv <- c(R0,rep(0,length(ages)-1))
  for (a in 1:(Amax-1)) {
    SSBv[a+1] <- SSBv[a]*exp(-M)
  }
  SSBv[Amax] <- SSBv[Amax-1]*exp(-M)/(1-exp(-M))
  SSB0 <- sum(SSBv*Wa*Mata*0.5)
  
  #### recruitment  ####
  if(rec==0){ # flat
    # 1st row = recruitment
    N[1,1] <- R0*exp(RecDevs[1] - 0.5*sigmaRec^2)
    for (t in 2:length(N[1,])) {
      N[1,t] <- R0*exp(RecDevs[t]-(0.5*sigmaRec^2))
    }
    
    for (a in 1:(Amax-1)) {
      N[a+1,1] <- N[a,1]*exp(-Zmat[a,1])
    }
    
    for (t in 1:(length(years)-1)) {
      for (a in 1:(Amax-1)) {
        N[a+1,t+1] <- N[a,t]*exp(-Zmat[a,t])
      }
      N[Amax,t+1] <- N[Amax-1,t]*exp(-Zmat[Amax-1,t]) + N[Amax, t] * exp(-Zmat[Amax, t])
    }
  } else if(rec==1){ # Beverton and Holt

    # 1st row = recruitment
    N[1,1] <- R0*exp(RecDevs[1] - 0.5*sigmaRec^2)
    for (a in 1:(Amax-1)) {
      N[a+1,1] <- N[a,1]*exp(-Zmat[a,1])
    }
    
    for (t in 1:(length(years)-1)) {
      ssbly <- sum(N[,t]*Wa*Mata*0.5)
      N[1,t+1] <- (4*h*R0*ssbly)/(SSB0*(1-h)+ssbly*(5*h-1))*exp(RecDevs[t+1] - 0.5*sigmaRec^2)
      
      for (a in 1:(Amax-1)) {
        N[a+1,t+1] <- N[a,t]*exp(-Zmat[a,t])
      }
      N[Amax, t+1] <- N[Amax-1, t] * exp(-Zmat[Amax-1, t]) + N[Amax, t] * exp(-Zmat[Amax, t])
    }
    
  } else if(rec==2){ # Ricker's 
    
    # 1st row = recruitment
    N[1,1] <- R0*exp(RecDevs[1] - 0.5*sigmaRec^2)
    for (a in 1:(Amax-1)) {
      N[a+1,1] <- N[a,1]*exp(-Zmat[a,1])
    }
    
    for (t in 1:(length(years)-1)) {
      ssbly <- sum(N[,t]*Wa*Mata*0.5)
      N[1,t+1] <- R0*(ssbly/SSB0)*exp((log(5*h)/0.8)*(1-ssbly/SSB0))* exp(RecDevs[t+1] - 0.5*sigmaRec^2)
      
      for (a in 1:(Amax-1)) {
        N[a+1,t+1] <- N[a,t]*exp(-Zmat[a,t])
      }
      N[Amax, t+1] <- N[Amax-1, t] * exp(-Zmat[Amax-1, t]) + N[Amax, t] * exp(-Zmat[Amax, t])
    }
  }
  
  ### nll RecDevs ####
  nll.recdevs.i <- -(dnorm(RecDevs,0,sigmaRec,log=TRUE))
  nll.recdevs <- sum(nll.recdevs.i)
  
  ### predicted catch ####
  logPredC <- log(colSums(Fmat/Zmat*(1-exp(-Fmat-M))*N*Wa))
  
  ### nll catch ####
  nll.catch.i <- -(dnorm(log(catch), logPredC, sigmaC, log=TRUE))
  nll.catch <- sum(nll.catch.i)
  
  ### abundance index ####
  Ba <- N*Wa*Sa
  IndexPred <- colSums(q*Ba)
  IndexPred <- IndexPred[index.idx]
  ### nll index ####
  nll.index.i <- -(dnorm(log(index.obs), log(IndexPred), sigmaIdx, log=TRUE))
  nll.index <- sum(nll.index.i)
  
  ### lengths ####
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
  
  length.prop <- length.obs/rowSums(length.obs)
  length.ESS <- length.nsamps*5
  length.adjusted <- length.prop*length.ESS
  
  nll.length.i <- rep(0,nrow(Pmat))
  
  ### nll lengths ####
  if(ESS==FALSE){
    for (i in 1:length(nll.length.i)) {
      nll.length.i[i] <- -sum(dmultinom(x=as.vector(length.obs[i,]), prob=Pmat[i,], log=TRUE))
    }
    nll.length <- sum(nll.length.i)
  } else {
    for (i in 1:length(nll.length.i)) {
      nll.length.i[i] <- -sum(dmultinom(x=as.vector(length.adjusted[i,]), prob=Pmat[i,], log=TRUE))
    }
    nll.length <- sum(nll.length.i)
  }
  ### nll total ####
  ret <- sum(nll.recdevs,nll.catch,nll.index,nll.length)

  ### reports ####
  
  SSB <- colSums(N*Wa*Mata*0.5)
  Rec<- N[1, ]
  B<- colSums(N*Wa)
  F <- apply(Fmat, 2, max) 
  Bratio <- SSB/SSB0
  
  ADREPORT(SSB0)
  ADREPORT(SSB)
  ADREPORT(Bratio)
  ADREPORT(Rec)
  ADREPORT(B)
  ADREPORT(F)
  ADREPORT(R0)
  ADREPORT(Wa)
  ADREPORT(Mata)
  ADREPORT(Sa)
  ADREPORT(nll.recdevs)
  ADREPORT(nll.catch)
  ADREPORT(nll.index)
  ADREPORT(nll.length)
  ADREPORT(ret)
  ADREPORT(nll.catch.i)
  ADREPORT(nll.index.i)
  ADREPORT(nll.length.i)
  ADREPORT(logPredC)
  ADREPORT(IndexPred)
  ADREPORT(Pmat)
  #####  
  return(ret)
}

#### model fitting and quick diags ####
{starttime <- Sys.time()
mapping <- list(M=factor(NA), h=factor(NA))
obj <- MakeADFun(nll,par, random = "RecDevs", map = mapping)
fit <- nlminb(obj$par,obj$fn,obj$gr)
sdrep <- sdreport(obj)
endtime <- Sys.time()
runtime <- endtime-starttime
round(runtime,2)}

fit$par
fit$convergence
max(abs(obj$gr(fit$par)))
sdrep$pdHess

pl <- as.list(sdrep, "Est")
plr <- as.list(sdrep, "Est", report=TRUE)
plr_se <- as.list(sdrep, "Std. Error", report = TRUE)

#### jitter analysis ####
# adapting do_jitter() function from SPoRC package
do_jitter1 <- function (data, parameters, mapping, random = NULL, sd, n_jitter, 
          n_newton_loops, do_par, n_cores, par_vec = NULL) 
{
  jitter_all <- data.frame()
  obj <- RTMB::MakeADFun(func = nll, parameters = parameters, 
                         map = mapping, random = random, silent = TRUE)
  if (do_par == FALSE) {
    for (i in 1:n_jitter) {
      if (is.null(par_vec)) 
        jitter_pars <- obj$par + stats::rnorm(length(obj$par), 
                                              0, sd)
      else jitter_pars <- par_vec + stats::rnorm(length(obj$par), 
                                                 0, sd)
      optim <- stats::nlminb(jitter_pars, obj$fn, obj$gr, 
                             control = list(iter.max = 1e+05, eval.max = 1e+05, 
                                            rel.tol = 1e-15))
      try_improve <- tryCatch(expr = for (j in 1:n_newton_loops) {
        g = as.numeric(obj$gr(optim$par))
        h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
        optim$par = optim$par - solve(h, g)
        optim$objective = obj$fn(optim$par)
      }, error = function(e) {
        e
      }, warning = function(w) {
        w
      })
      obj$rep <- obj$report(obj$env$last.par.best)
      obj$sd_rep <- RTMB::sdreport(obj)
      plr1 <- as.list(obj$sd_rep, "Est", report=TRUE)
      jitter_ts_df <- data.frame(
        Year = data$years,
        SSB = as.numeric(plr1$SSB),
        Rec = as.numeric(plr1$Rec))
      jitter_ts_df$jitter  <-  i
      jitter_ts_df$Hessian  <-  obj$sd_rep$pdHess
      jitter_ts_df$jnLL  <-  optim$objective
      jitter_ts_df$Max_Gradient <-  max(abs(obj$sd_rep$gradient.fixed))
      
      jitter_all <- rbind(jitter_all, jitter_ts_df)
    }
  }
  if (do_par == TRUE) {
    future::plan(future::multisession, workers = n_cores)
    progressr::with_progress({
      p <- progressr::progressor(along = 1:n_jitter)
      jitter_all <- future.apply::future_lapply(1:n_jitter, 
                                                function(i) {
                                                  obj <- RTMB::MakeADFun(func = nll, 
                                                                         parameters = parameters, map = mapping, random = random, 
                                                                         silent = TRUE)
                                                  if (is.null(par_vec)) 
                                                    jitter_pars <- obj$par + stats::rnorm(length(obj$par), 
                                                                                          0, sd)
                                                  else jitter_pars <- par_vec + stats::rnorm(length(obj$par), 
                                                                                             0, sd)
                                                  optim <- stats::nlminb(jitter_pars, obj$fn, 
                                                                         obj$gr, control = list(iter.max = 1e+05, 
                                                                                                eval.max = 1e+05, rel.tol = 1e-15))
                                                  try_improve <- tryCatch({
                                                    for (j in 1:n_newton_loops) {
                                                      g <- as.numeric(obj$gr(optim$par))
                                                      h <- optimHess(optim$par, fn = obj$fn, 
                                                                     gr = obj$gr)
                                                      optim$par <- optim$par - solve(h, g)
                                                      optim$objective <- obj$fn(optim$par)
                                                    }
                                                  }, error = function(e) e, warning = function(w) w)
                                                  obj$rep <- obj$report(obj$env$last.par.best)
                                                  obj$sd_rep <- RTMB::sdreport(obj)
                                                  plr1 <- as.list(obj$sd_rep, "Est", report=TRUE)
                                                  jitter_ts_df <- data.frame(
                                                    Year = data$years,
                                                    SSB = as.numeric(plr1$SSB),
                                                    Rec = as.numeric(plr1$Rec))
                                                  jitter_ts_df$jitter  <-  i
                                                  jitter_ts_df$Hessian  <-  obj$sd_rep$pdHess
                                                  jitter_ts_df$jnLL  <-  optim$objective
                                                  jitter_ts_df$Max_Gradient <-  max(abs(obj$sd_rep$gradient.fixed))
                                                  p()
                                                  jitter_ts_df
                                                }, future.seed = TRUE) %>% bind_rows()
      future::plan(future::sequential)
    })
  }
  return(jitter_all)
}

jit <- do_jitter1(data = dat, parameters = par, mapping = mapping, random = "RecDevs", 
           sd = 0.2, n_jitter = 50, n_newton_loops = 1, do_par = TRUE, n_cores = 6)
str(jit)
write.table(jit,"outputs/jitter_results.txt", row.names = FALSE, sep = "\t", col.names = TRUE)
table(round(jit$jnLL, 3))
sum(jit$Hessian == FALSE)
range(jit$jnLL)

#### likelihood profiles (M and h) ####
#### retrospective ####
#### hindcast ####
#### plots ####

# 3 panel catch, index fit and lengths fit
{
  #x11()
  #pdf("outputs/model_fits.pdf", width = 8, height = 11)
  wd <- 3000
  ht <- 4000
  re <- 300
  png("outputs/model_fits.png", width = wd, height = ht, res = re)
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  # catch timeseries
  plot(dat$years, dat$catch, pch = 20, type = "h",lwd = 10,
       xlab = "Year", ylab = "Catch (Metric tons)", col = "gray70",
       main = "Catch time series")
  # cpue residuals
  plot(dat$index.yrs, dat$index.obs, pch = 20,
       xlab = "Year", ylab = "Abundance index",
       main = "Observed vs Predicted abundance index", ylim = c(0, 12))
  lines(dat$index.yrs, plr$IndexPred, col = "red", lwd = 2)
  legend("topright", legend = c("Observed", "Predicted"), 
         col = c("black", "red"), pch = c(20, NA), lty = c(NA, 1), lwd = c(NA, 2))
  
  # length residuals
  plot(dat$length.bins, 
       colSums(dat$length.obs)/sum(colSums(dat$length.obs)), 
       type = "h", lwd = 6, col = "gray70",
       xlab = "Length (cm)", ylab = "Proportion", 
       main = "Aggregated Length Composition",
       ylim = c(0, 0.05))
  lines(dat$length.bins, 
        colSums(matrix(plr$Pmat, nrow = length(dat$length.idx),ncol = length(dat$length.bins))
                                 *dat$length.nsamps)/sum(dat$length.nsamps)
        , col = "red", lwd = 2.5, type = "b", pch = 16, cex = 0.6)
  legend("topright", legend = c("Observed", "Predicted"), 
         col = c("gray70", "red"), lwd = c(6, 2.5), lty = c(1, 1), pch = c(NA, 16))
  graphics.off()
}

# 2 panel catch residuals and index residuals
{ 
  #x11()
  #pdf("outputs/residuals_plots.pdf", width = 8, height = 8)
  wd <- 3000
  ht <- 4000
  re <- 300
  png("outputs/residuals_plots.png", width = wd, height = ht, res = re)
  par(mfrow=c(2,1))
  
  # catch residuals
  resid.catch <- log(dat$catch) - plr$logPredC
  plot(dat$years, resid.catch, pch = 20,
       xlab = "Year", ylab = "Log Residuals",
       main = "Catch Residuals", ylim = c(-1, 1))
  abline(h = 0, col = "red", lty = 1)
  
  # cpue residuals
  resid.index <- log(dat$index.obs) - log(plr$IndexPred)
  plot(dat$index.yrs, resid.index, pch=20, ylim=c(-1,1), main="CPUE Residuals",
       xlab = "Year", ylab = "Log Residuals",)
  abline(h=0, col="red", lty=1)
  graphics.off()
}

# 4 panel biomass, bratio, RecDevs and F
{
  wd <- 4000
  ht <- 3000
  re <- 300
  #x11()
  png("outputs/results_2x2.png", width = wd, height = ht, res = 300)
  par(mfrow=c(2,2))
  plot(dat$years,plr$B, ylim = c(0,max(plr$B + 1.96 * plr_se$B)), pch=20,
       ylab = "Biomass (t)", main = "Total Biomass")
  lines(dat$years,plr$B)
  lines(dat$years,plr$B - 1.96 * plr_se$B, lty=2)
  lines(dat$years,plr$B + 1.96 * plr_se$B, lty=2)
  
  plot(dat$years,plr$Bratio, 
       ylim = c(0,as.numeric(number(max(plr$Bratio + 1.96 * plr_se$Bratio),0.1))), pch=20,
       xlab = "Year", ylab = expression(Depletion ~ (SSB / SSB[0])), main = "Stock Depletion (Bratio)")
  lines(dat$years,plr$Bratio)
  lines(dat$years,plr$Bratio - 1.96 * plr_se$Bratio, lty=2)
  lines(dat$years,plr$Bratio + 1.96 * plr_se$Bratio, lty=2)
  
  plot(dat$years,plr$Rec/1e6, pch=20,
       xlab = "Year", ylab = "Age-0 recruits (Millions)", main = "Recruitment")
  lines(dat$years,plr$Rec/1e6)
  plot(dat$years,plr$F, pch=20, ylim = c(0,max(plr$F + 1.96*plr_se$F)),
       xlab = "Year", ylab = "F", main = "Fishing Mortality")
  arrows(dat$years, plr$F - 1.96*plr_se$F,
         dat$years, plr$F + 1.96*plr_se$F,
         angle=0, code=3, length=0.02)
  graphics.off()

}

# Jitter analysis
{
  wd <- 4000
  ht <- 3000
  re <- 300
  png("outputs/jitter.png", width = wd, height = ht, res = 300)
  plot(jit$jitter, jit$jnLL, pch=20, col=ifelse(jit$jnLL > (min(jit$jnLL*1.05)), "darkred", "darkblue"),
       xlab = "Jitter run", ylab = "Negative log likelihood value",
       ylim = c(min(jit$jnLL) - 0.01, median(jit$jnLL) + 0.01),
       main = "Jitter NLL Distribution"); abline(h=mean(jit$jnLL))
  legend("topright", legend = c("Converged (Optimal)", "Local minimum (Outlier)", "Mean NLL"), 
         col = c("darkblue", "darkred", "black"), 
         pch = c(20, 20, NA), lty = c(NA, NA, 1), bty = "n")
  graphics.off()
}



















