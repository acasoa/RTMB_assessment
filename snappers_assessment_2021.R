setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(RTMB)
library(SPoRC)
library(scales)
library(dplyr)
library(parallel)
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
  ESS=FALSE,
  n_fish_fleets=1
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

#### adapted functions from SPoRC ####
## RPS
Get_Reference_Points1 <- function (data, rep, SPR_x = NULL, t_spwn = 0,sex_ratio_f = 0.5, calc_rec_st_yr = 1, 
                                   rec_age = 1, type = "single_region", what="BH_MSY", n_avg_yrs = 1, 
                                   local_bh_msy_newton_steps = 6, M=NULL, h_trans=NULL) 
{
  f_ref_pt <- vector()
  b_ref_pt <- vector()
  virgin_b_ref_pt <- vector()
  n_yrs <- length(data$years)
  avg_yrs <- (n_yrs - n_avg_yrs + 1):n_yrs
  if (type == "single_region") {
    if (!what %in% c("SPR", "BH_MSY")) 
      stop("what is not correctly specified! Should be SPR, BH_MSY for type = single_region")
    data_list <- list()
    n_ages <- length(data$ages)
    n_years <- length(data$years)
    data_list$t_spwn <- t_spwn
    data_list$F_fract_flt <- 1#rep$F[1, n_years, ]/sum(rep$F[1, 
    #                                    n_years, ])
    fish_sel_avg <- rep$Sa#apply(rep$fish_sel[, avg_yrs, , 1, , 
    #            drop = FALSE], c(1, 3, 4, 5), mean)
    data_list$fish_sel <- array(fish_sel_avg, dim = c(n_ages, 
                                                      data$n_fish_fleets))
    natmort_avg <- #apply(rep$natmort[, avg_yrs, , 1, drop = FALSE], 
      #c(1, 3, 4), mean)
      data_list$natmort <- as.vector(M)
    #WAA_avg <- apply(data$WAA[, avg_yrs, , 1, drop = FALSE], 
    #                 c(1, 3, 4), mean)
    data_list$WAA <- rep$Wa
    #MatAA_avg <- apply(data$MatAA[, avg_yrs, , 1, drop = FALSE], 
    #     c(1, 3, 4), mean)
    data_list$MatAA <- rep$Mata#MatAA_avg
    data_list$sex_ratio_f <- sex_ratio_f
    
    single_region_SPR1 <- function (par_list){
      
      "c" <- RTMB::ADoverload("c")
      "[<-" <- RTMB::ADoverload("[<-")
      RTMB::getAll(par_list, data_list)
      n_ages = dim(fish_sel)[1]
      F_x = exp(log_F_x)
      SB_age = Nspr = array(0, dim = c(2, n_ages))
      Nspr[, 1] = 1 * sex_ratio_f
      for (j in 2:(n_ages - 1)) {
        Nspr[1, j] = Nspr[1, j - 1] * exp(-1 * M)
        Nspr[2, j] = Nspr[2, j - 1] * exp(-1 * (M + 
                                                  sum(F_fract_flt * F_x * fish_sel[j - 1, ])))
      }
      Nspr[1, n_ages] = Nspr[1, n_ages - 1] * exp(-1 * M)/(1 - exp(-1 * M))
      Nspr[2, n_ages] = Nspr[2, n_ages - 1] * exp(-1 * M + sum(F_fract_flt * F_x * fish_sel[n_ages - 1, ]))/(1 -exp(-1 * (M + sum(F_fract_flt * F_x * fish_sel[n_ages, ]))))
      for (j in 1:n_ages) {
        SB_age[1, j] = Nspr[1, j] * WAA[j] * MatAA[j] * exp(-t_spwn * 
                                                              M)
        SB_age[2, j] = Nspr[2, j] * WAA[j] * MatAA[j] * exp(-t_spwn * 
                                                              (M + sum(F_fract_flt * F_x * fish_sel[j, 
                                                              ])))
      }
      SB0 = sum(SB_age[1, ])
      SB_F_x = sum(SB_age[2, ])
      SPR = SB_F_x/SB0
      sprpen = 100 * (SPR - SPR_x)^2
      RTMB::REPORT(SB_age)
      RTMB::REPORT(Nspr)
      RTMB::REPORT(SB0)
      RTMB::REPORT(SB_F_x)
      RTMB::REPORT(F_x)
      return(sprpen)
    }
    single_region_BH_Fmsy1 <- function (par_list) {
      "c" <- RTMB::ADoverload("c")
      "[<-" <- RTMB::ADoverload("[<-")
      RTMB::getAll(par_list, data_list)
      n_ages = dim(fish_sel)[1]
      Fmsy = exp(log_Fmsy)
      SB_age = Nspr = array(0, dim = c(2, n_ages))
      CAA = rep(0, n_ages)
      Nspr[, 1] = 1 * sex_ratio_f
      for (j in 2:(n_ages - 1)) {
        Nspr[1, j] = Nspr[1, j - 1] * exp(-1 * M)
        Nspr[2, j] = Nspr[2, j - 1] * exp(-1 * (M + 
                                                  sum(F_fract_flt * Fmsy * fish_sel[j - 1, ])))
      }
      Nspr[1, n_ages] = Nspr[1, n_ages - 1] * exp(-1 * M)/(1 - exp(-1 * M))
      Nspr[2, n_ages] = Nspr[2, n_ages - 1] * (exp(-1 * (M + sum(F_fract_flt * Fmsy * fish_sel[n_ages - 1, ])))/
                                                 (1 - exp(-1 * (M + sum(F_fract_flt * Fmsy * 
                                                                          fish_sel[n_ages, ])))))
      for (j in 1:n_ages) {
        tmp_F = sum(F_fract_flt * Fmsy * fish_sel[j, ])
        tmp_Z = tmp_F + M
        SB_age[1, j] = Nspr[1, j] * WAA[j] * MatAA[j] * exp(-t_spwn * 
                                                              M)
        SB_age[2, j] = Nspr[2, j] * WAA[j] * MatAA[j] * exp(-t_spwn * 
                                                              tmp_Z)
        CAA[j] = Nspr[2, j] * (tmp_F/tmp_Z) * (1 - exp(-tmp_Z))
      }
      SBPR_0 = sum(SB_age[1, ])
      SBPR_F = sum(SB_age[2, ])
      Req = R0 * ((4 * h * SBPR_F) - (1 - h) * SBPR_0)/((5 * h - 
                                                           1) * SBPR_F)
      Yield = sum(CAA * WAA) * Req
      Bmsy = SBPR_F * Req
      B0 = SBPR_0 * R0
      obj_fun = -Yield
      RTMB::REPORT(SB_age)
      RTMB::REPORT(Nspr)
      RTMB::REPORT(SBPR_0)
      RTMB::REPORT(SBPR_F)
      RTMB::REPORT(Fmsy)
      RTMB::REPORT(Yield)
      RTMB::REPORT(Bmsy)
      RTMB::REPORT(B0)
      RTMB::REPORT(Req)
      return(obj_fun)
    }
    
    if (what == "SPR") {
      data_list$SPR_x <- SPR_x
      par_list <- list()
      par_list$log_F_x <- log(0.1)
      
      obj <- RTMB::MakeADFun(func=single_region_SPR1, 
                             parameters = par_list, map = NULL, silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr)
      obj$rep <- obj$report(obj$env$last.par.best)
      f_ref_pt[1] <- obj$rep$F_x
      b_ref_pt[1] <- obj$rep$SB_F_x * rep$R0
      virgin_b_ref_pt[1] <- obj$rep$SB0 * rep$R0
    }
    if (what == "BH_MSY") {
      data_list$h <- h_trans
      data_list$R0 <- rep$R0
      par_list <- list()
      par_list$log_Fmsy <- log(0.1)
      obj <- RTMB::MakeADFun(single_region_BH_Fmsy1, 
                             parameters = par_list, map = NULL, 
                             silent = TRUE)
      obj$optim <- stats::nlminb(obj$par, obj$fn, obj$gr)
      obj$rep <- obj$report(obj$env$last.par.best)
      f_ref_pt[1] <- obj$rep$Fmsy
      b_ref_pt[1] <- obj$rep$Bmsy
      virgin_b_ref_pt[1] <- obj$rep$B0
    }
  }
  return(list(f_ref_pt = f_ref_pt, b_ref_pt = b_ref_pt, 
              virgin_b_ref_pt = virgin_b_ref_pt))
}

## Jitter
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

## likelihood profiles
do_likelihood_profile1 <- function (data, parameters, mapping, random = NULL, what, idx = NULL, 
                                    min_val, max_val, inc = 0.05, do_par = TRUE, n_cores = NULL) 
{
  if (min_val > max_val) {
    stop("`min_val` is greater than `max_val`. This likely occurred because you are profiling a log-transformed parameter. 
         Try swapping the values: use the current `min_val` as `max_val`, and vice versa.")
  }
  vals <- seq(min_val, max_val, inc)
  jnLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, 
                                                                    NULL))
  rec_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, 
                                                                       NULL))
  Catch_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, 
                                                                         NULL))
  FishLen_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, 
                                                                           NULL))
  FishIdx_nLL <- matrix(NA, nrow = length(vals), ncol = 1, dimnames = list(vals, 
                                                                           NULL))
  if (do_par == FALSE) {
    for (j in 1:length(vals)) {
      if (!is.null(dim(parameters[[what]]))) {
        parameters[[what]] <- do.call(`[<-`, c(list(parameters[[what]]), 
                                               idx, list(vals[j])))
        counter <- 1
        map_parameter <- do.call(`[<-`, c(list(parameters[[what]]), 
                                          idx, list(NA)))
        non_na <- which(!is.na(map_parameter))
        for (i in 1:length(non_na)) {
          map_parameter[non_na[i]] <- counter
          counter <- counter + 1
        }
        mapping[[what]] <- factor(map_parameter)
      }else {
        parameters[[what]] <- vals[j]
        mapping[[what]] <- factor(NA)
      }
      SPoRC_rtmb_model <- RTMB::MakeADFun(func = nll 
                                          , parameters = parameters, map = mapping, 
                                          random = random, silent = TRUE)
      tryCatch({
        SPoRC_optim <- stats::nlminb(SPoRC_rtmb_model$par, 
                                     SPoRC_rtmb_model$fn, SPoRC_rtmb_model$gr)#, control = list(iter.max = 1e+06, 
        #               eval.max = 1e+06, rel.tol = 1e-15))
        #report <- SPoRC_rtmb_model$report(SPoRC_rtmb_model$env$last.par.best)
        SPoRC_rtmb_model$sd_rep <- RTMB::sdreport(SPoRC_rtmb_model)
        plr1 <- as.list(SPoRC_rtmb_model$sd_rep, "Est", report=TRUE)
        jnLL[j, 1] <- plr1$ret
        rec_nLL[j, 1] <- plr1$nll.recdevs
        Catch_nLL[j, 1] <- plr1$nll.catch
        FishLen_nLL[j, 1] <- plr1$nll.length
        FishIdx_nLL[j, 1] <- plr1$nll.index
        
        print(paste("Likelihood profile is at:", round(j/length(vals) * 
                                                         100, 2), "%"))
      }, error = function(e) {
        message("Failed to optimize: ", e$message)
      })
    }
  }
  if (do_par == TRUE) {
    future::plan(future::multisession, workers = n_cores)
    progressr::with_progress({
      p <- progressr::progressor(along = 1:length(vals))
      profile_results <- future.apply::future_lapply(1:length(vals), 
                                                     function(j) {
                                                       local_parameters <- parameters
                                                       local_mapping <- mapping
                                                       result <- list(j = j, prof_val = vals[j], success = FALSE, 
                                                                      jnLL = NA, rec_nLL = NA, Catch_nLL = data.frame(),
                                                                      FishLen_nLL = data.frame(), FishIdx_nLL = data.frame())
                                                       if (!is.null(dim(local_parameters[[what]]))) {
                                                         local_parameters[[what]] <- do.call(`[<-`, 
                                                                                             c(list(local_parameters[[what]]), idx, 
                                                                                               list(vals[j])))
                                                         counter <- 1
                                                         map_parameter <- do.call(`[<-`, c(list(local_parameters[[what]]), 
                                                                                           idx, list(NA)))
                                                         non_na <- which(!is.na(map_parameter))
                                                         for (i in 1:length(non_na)) {
                                                           map_parameter[non_na[i]] <- counter
                                                           counter <- counter + 1
                                                         }
                                                         local_mapping[[what]] <- factor(map_parameter)
                                                       } else {
                                                         local_parameters[[what]] <- vals[j]
                                                         local_mapping[[what]] <- factor(NA)
                                                       }
                                                       tryCatch({
                                                         SPoRC_rtmb_model <- RTMB::MakeADFun(func = nll 
                                                                                             , parameters = local_parameters, 
                                                                                             map = local_mapping, 
                                                                                             random = random, silent = TRUE)
                                                         SPoRC_optim <- stats::nlminb(SPoRC_rtmb_model$par, 
                                                                                      SPoRC_rtmb_model$fn, SPoRC_rtmb_model$gr)
                                                         #report <- SPoRC_rtmb_model$report(SPoRC_rtmb_model$env$last.par.best)
                                                         SPoRC_rtmb_model$sd_rep <- RTMB::sdreport(SPoRC_rtmb_model)
                                                         plr1 <- as.list(SPoRC_rtmb_model$sd_rep, "Est", report=TRUE)
                                                         result$jnLL <- plr1$ret
                                                         result$rec_nLL <- plr1$nll.recdevs
                                                         result$Catch_nLL <- plr1$nll.catch
                                                         result$FishLen_nLL <- plr1$nll.length
                                                         result$FishIdx_nLL <- plr1$nll.index
                                                         result$success <- TRUE
                                                       }, error = function(e) {
                                                         result$error_msg <- e$message
                                                       })
                                                       p(sprintf("j=%g", j))
                                                       return(result)
                                                     }, future.seed = TRUE)
    })
    for (result in profile_results) {
      if (result$success) {
        j <- result$j
        jnLL[j, 1] <- result$jnLL
        rec_nLL[j, 1] <- result$rec_nLL
        Catch_nLL[j, 1] <- result$Catch_nLL
        FishLen_nLL[j, 1] <- result$FishLen_nLL
        FishIdx_nLL[j, 1] <- result$FishIdx_nLL
      }
      else {
        message("Failed iteration ", result$j, " (value = ", 
                result$prof_val, "): ", result$error_msg)
      }
    }
    future::plan(future::sequential)
    message("Parallel processing completed.")
  }
  jnLL_df <- reshape2::melt(jnLL) %>% dplyr::select(-Var2) %>% 
    dplyr::rename(prof_val = Var1) %>% dplyr::mutate(type = "jnLL")
  rec_nLL_df <- reshape2::melt(rec_nLL) %>% dplyr::select(-Var2) %>% 
    dplyr::rename(prof_val = Var1) %>% dplyr::mutate(type = "RecPen")
  Catch_nLL_df <- reshape2::melt(Catch_nLL) %>% dplyr::select(-Var2) %>% 
    dplyr::rename(prof_val = Var1) %>% dplyr::mutate(type = "Catch")
  FishLen_nLL_df <- reshape2::melt(FishLen_nLL) %>% dplyr::select(-Var2) %>% 
    dplyr::rename(prof_val = Var1) %>% dplyr::mutate(type = "FishLen")
  FishIdx_nLL_df <- reshape2::melt(FishIdx_nLL) %>% dplyr::select(-Var2) %>% 
    dplyr::rename(prof_val = Var1) %>% dplyr::mutate(type = "FishIdx")
  
  profile_list <- list(jnLL_df = jnLL_df, rec_nLL_df = rec_nLL_df, 
                       Catch_nLL_df = Catch_nLL_df, FishLen_nLL_df = FishLen_nLL_df, 
                       FishIdx_nLL_df = FishIdx_nLL_df)
  return(profile_list)
}

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
random <- "RecDevs"
obj <- MakeADFun(nll,par, random = random, map = mapping, silent = TRUE)
fit <- nlminb(obj$par,obj$fn,obj$gr)
sdrep <- sdreport(obj, silent = TRUE)
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

BH_RP <- Get_Reference_Points1(data = dat,rep = plr,SPR_x = 0.4, M = par$M, h_trans = par$h)
SPR_RP <- Get_Reference_Points1(data = dat,rep = plr,SPR_x = 0.4, M = par$M, h_trans = par$h, what = "SPR")

saveRDS(list(BH_RP = BH_RP, SPR_RP = SPR_RP, fit=fit$par, 
             plr= plr[c("SSB", "Rec", "B", "F", "Bratio", "SSB0", "R0")],
             plr_se= plr_se[c("SSB", "Rec", "B", "F", "Bratio", "SSB0", "R0")]), "outputs/outputs.rds")

#### jitter analysis ####
# adapting do_jitter() function from SPoRC package

n_cores <- 6; sd <- 0.2; n_jitter <- 50
jit <- do_jitter1(data = dat, parameters = par, mapping = mapping, random = random, 
           sd = sd, n_jitter = n_jitter, n_newton_loops = 1, do_par = TRUE, n_cores = n_cores)
str(jit)
table(round(jit$jnLL, 3))
sum(jit$Hessian == FALSE)
range(jit$jnLL)
saveRDS(jit, "outputs/jitter_results.rds")

#### likelihood profiles (M and h) ####
# adapting do_likelihood_profile() function from SPoRC package
what <- c("M","h"); min_val <- c(0.15,0.5); max_val <- c(0.25,0.9); inc <- c(0.01,0.05); n_cores <- 6
profiles <- list()

for (i in 1:length(what)) {
  cat("Profiling", what[i], "\n")

  profiles[[what[i]]] <- do_likelihood_profile1(data = dat,
                                                parameters = par,
                                                mapping = mapping,
                                                random = random,
                                                what = what[i],
                                                min_val = min_val[i],
                                                max_val = max_val[i],
                                                inc = inc[i],
                                                do_par = TRUE,
                                                n_cores = n_cores
  )
}

saveRDS(profiles,"outputs/likelihood_profiles.rds")

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
       xlab = "Year",ylab = "Biomass (t)", main = "Total Biomass")
  lines(dat$years,plr$B)
  lines(dat$years,plr$B - 1.96 * plr_se$B, lty=2)
  lines(dat$years,plr$B + 1.96 * plr_se$B, lty=2)
  
  plot(dat$years,plr$Bratio, 
       ylim = c(0,as.numeric(number(max(plr$Bratio + 1.96 * plr_se$Bratio),0.1))), pch=20,
       xlab = "Year", ylab = expression(Depletion ~ (SSB / SSB[0])), main = "Stock Depletion (Bratio)")
  lines(dat$years,plr$Bratio)
  lines(dat$years,plr$Bratio - 1.96 * plr_se$Bratio, lty=2)
  lines(dat$years,plr$Bratio + 1.96 * plr_se$Bratio, lty=2)
  abline(h=BH_RP$b_ref_pt/BH_RP$virgin_b_ref_pt, col="red", lty=2)
  text(x = min(dat$years) + 4,y = BH_RP$b_ref_pt/BH_RP$virgin_b_ref_pt,
       labels = expression("Management Target"), pos = 3,col = "black",cex = 0.7)
  
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

# likelihood profiles
{
  if(length(profiles)>=2){
    wd <- 4000
    ht <- 3000
    re <- 300
    for (i in 1:length(profiles)) {
      png(paste0("outputs/likelihood_profiles_",
                 ifelse(i==1,"NatM","Steepbh"),".png"), width = wd, height = ht, res = 300)
      
      plot(as.numeric(seq(min_val[[i]], max_val[[i]], inc[[i]])),
           profiles[[i]]$jnLL_df$value - min(profiles[[i]]$jnLL_df$value), 
           type = "b", pch = 19, lwd = 2, col = "black", cex = 0.8,
           ylim = c(0, max(profiles[[i]]$jnLL_df$value - min(profiles[[i]]$jnLL_df$value),
                           profiles[[i]]$Catch_nLL_df$value - min(profiles[[i]]$Catch_nLL_df$value),
                           profiles[[i]]$FishLen_nLL_df$value - min(profiles[[i]]$FishLen_nLL_df$value),
                           profiles[[i]]$FishIdx_nLL_df$value - min(profiles[[i]]$FishIdx_nLL_df$value),
                           profiles[[i]]$rec_nLL_df$value - min(profiles[[i]]$rec_nLL_df$value)) * 1.1),
           xlab = ifelse(i==1,"Natural Mortality (M)","Steepness (h)"),
           ylab = expression(Delta ~ "Negative Log-Likelihood"),
           main = "Likelihood Profile Components")
      
      lines(seq(min_val[[i]], max_val[[i]], inc[[i]]), profiles[[i]]$Catch_nLL_df$value - min(profiles[[i]]$Catch_nLL_df$value), type = "b", pch = 16, col = "red", cex = 0.8, lty = 2)
      lines(seq(min_val[[i]], max_val[[i]], inc[[i]]), profiles[[i]]$FishLen_nLL_df$value - min(profiles[[i]]$FishLen_nLL_df$value),   type = "b", pch = 17, col = "blue", cex = 0.8, lty = 2)
      lines(seq(min_val[[i]], max_val[[i]], inc[[i]]), profiles[[i]]$FishIdx_nLL_df$value - min(profiles[[i]]$FishIdx_nLL_df$value),   type = "b", pch = 18, col = "green4", cex = 0.8, lty = 2)
      lines(seq(min_val[[i]], max_val[[i]], inc[[i]]), profiles[[i]]$rec_nLL_df$value - min(profiles[[i]]$rec_nLL_df$value),   type = "b", pch = 15, col = "orange", cex = 0.8, lty = 2)
  
      legend("topright", legend = c("Total", "Catch", "Length", "Index", "RecDevs"),
             col = c("black", "red", "blue", "green4", "orange"),
             pch = c(19, 16, 17, 18, 15), lty = c(1, 2, 2, 2, 2), cex = 0.8, bty = "n")
    
      graphics.off()
    }
  }else{
    png(paste0("outputs/likelihood_profiles_",
               ifelse(names(profiles[1])=="M","NatM","Steepbh"),".png"), width = wd, height = ht, res = 300)
    ifelse(names(profiles[1])=="M",j <- 1,j <- 2); p=1
    plot(as.numeric(seq(min_val[[j]], max_val[[j]], inc[[j]])),
         profiles[[p]]$jnLL_df$value - min(profiles[[p]]$jnLL_df$value),
         type = "b", pch = 19, lwd = 2, col = "black", cex = 0.8,
         ylim = c(0, max(profiles[[p]]$jnLL_df$value - min(profiles[[p]]$jnLL_df$value),
                         profiles[[p]]$Catch_nLL_df$value - min(profiles[[p]]$Catch_nLL_df$value),
                         profiles[[p]]$FishLen_nLL_df$value - min(profiles[[p]]$FishLen_nLL_df$value),
                         profiles[[p]]$FishIdx_nLL_df$value - min(profiles[[p]]$FishIdx_nLL_df$value),
                         profiles[[p]]$rec_nLL_df$value - min(profiles[[p]]$rec_nLL_df$value)) * 1.1),
         xlab = ifelse(j==1,"Natural Mortality (M)","Steepness (h)"),
         ylab = expression(Delta ~ "Negative Log-Likelihood"),
         main = "Likelihood Profile Components")
    
    lines(seq(min_val[[j]], max_val[[j]], inc[[j]]), profiles[[p]]$Catch_nLL_df$value - min(profiles[[p]]$Catch_nLL_df$value), type = "b", pch = 16, col = "red", cex = 0.8, lty = 2)
    lines(seq(min_val[[j]], max_val[[j]], inc[[j]]), profiles[[p]]$FishLen_nLL_df$value - min(profiles[[p]]$FishLen_nLL_df$value),   type = "b", pch = 17, col = "blue", cex = 0.8, lty = 2)
    lines(seq(min_val[[j]], max_val[[j]], inc[[j]]), profiles[[p]]$FishIdx_nLL_df$value - min(profiles[[p]]$FishIdx_nLL_df$value),   type = "b", pch = 18, col = "green4", cex = 0.8, lty = 2)
    lines(seq(min_val[[j]], max_val[[j]], inc[[j]]), profiles[[p]]$rec_nLL_df$value - min(profiles[[p]]$rec_nLL_df$value),   type = "b", pch = 15, col = "orange", cex = 0.8, lty = 2)
    
    legend("topright", legend = c("Total", "Catch", "Length", "Index", "RecDevs"),
           col = c("black", "red", "blue", "green4", "orange"),
           pch = c(19, 16, 17, 18, 15), lty = c(1, 2, 2, 2, 2), cex = 0.8, bty = "n")
    
    graphics.off()
  }
}

#####