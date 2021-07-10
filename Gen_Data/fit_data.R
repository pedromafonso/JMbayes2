setwd( "C:/Users/pedro/Documents/GitHub/JMbayes2-RE")
ref <- "18"
desc <- "frailty sampler"
n_data <- 100L
n <- 500L
scale <- "gap"
################################################################################
## fit data
dataL <- readRDS(file = paste0("Gen_Data/", ref,"/dataL_", ref,".rds"))
n_data <- length(dataL)
n_cores <- max(parallel::detectCores() - 1L, 1L)
#
n_chains <- 3L
ncl_in <- min(n_chains, n_cores) # number of inner clusters (per outer cluster)
ncl_out <- max(floor(n_cores / ncl_in), 1) # number of outer clusters
outer <- function(i, ncl_in){
  tic2 <- Sys.time()
  # fit data
  lme_fit <- lme(y ~ ns(time, k =  c(1, 3), B = c(0, 7) ), 
                 random = list(id = pdDiag(form = ~ ns(time, k = c(1, 3), 
                                                       B = c(0, 7)))),
                 data = dataL[[i]]$long, 
                 control = lmeControl(opt = "optim", niterEM = 45))
  cox_fit <- coxph(Surv(tstart, tstop, status) ~ (group + age) * strata(strata),
                   data = dataL[[i]]$rec_strt)
  jm_fit <- jm(cox_fit, lme_fit, time_var = "time", 
               functional_forms = list("y" = ~ value(y) * strata),
               cores = ncl_in, recurrent = scale)
  timer <- difftime(Sys.time(), tic2, units = "mins")
  # get estimates
  lme_betas <- fixef(lme_fit)
  lme_vcov <- vcov(lme_fit)[lower.tri(vcov(lme_fit), diag = TRUE)]
  lme_sigma <- sigma(lme_fit)
  lme_est <- list(betas = lme_betas, vcov = lme_vcov, sigma = lme_sigma)
  jm_betas <- jm_fit$statistics$Mean$betas1
  jm_vcov  <- jm_fit$statistics$Mean$D
  jm_sigma <- jm_fit$statistics$Mean$sigmas
  jm_alpha <- jm_fit$statistics$Mean$alphas
  jm_gammas <- jm_fit$statistics$Mean$gammas
  jm_alphaF <- jm_fit$statistics$Mean$alphaF
  jm_sigmaF <- jm_fit$statistics$Mean$sigmaF
  jm_frailty <- jm_fit$statistics$Mean$frailty
  jm_est <- list(betas = jm_betas, vcov = jm_vcov, sigma = jm_sigma, 
                 alpha = jm_alpha, gammas = jm_gammas, alphaF = jm_alphaF,
                 sigmaF = jm_sigmaF, frailty = jm_frailty)
  # return
  list(lme_est = lme_est, jm_est = jm_est, timer = timer)
}
tic1 <- Sys.time()
cl_out <- parallel::makeCluster(ncl_out)
invisible(parallel::clusterEvalQ(cl_out, library(JMbayes2))) # load packages in each cl
parallel::clusterExport(cl_out, list("dataL", "scale")) # load vars in each cl
res <- parallel::parLapply(cl_out, seq_len(n_data), outer, ncl_in)
parallel::stopCluster(cl_out)
toc <- Sys.time()
saveRDS(res, file = paste0("Gen_Data/", ref, "/estimates_", ref,".rds"))
diff_time <- sapply(res, "[[", "timer")
source("Gen_Data/sim_funs.R") # to use par_diag()
par_diag(n_data, n_chains, ncl_in, ncl_out, n_cores, diff_time)
dur_min <- difftime(toc, tic1, units = "mins")
saveRDS(dur_min, file = paste0("Gen_Data/", ref, "/dur_min_", ref,".rds"))
invisible(file.copy(from = "Gen_Data/fit_data.R", 
                    to = paste0("Gen_Data/", ref, "/fit_data.R"), 
                    overwrite = TRUE))
RPushbullet::pbPost("note", 
                    title = paste0("RecData: fit complete (", 
                                   round(dur_min/60, 2), " h)"))
print(round(difftime(toc, tic1), 2))