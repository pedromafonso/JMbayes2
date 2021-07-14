# changed to save the full mcmc for the alphaF and sigmaF
setwd( "C:/Users/pedro/Documents/GitHub/JMbayes2-RE")
ref <- "37"
scale <- "gap"
n_iter <- 7000L
n_burnin <- 1000L
n_thin <- 1L
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
               cores = ncl_in, recurrent = scale,
               n_iter = n_iter, n_burnin = n_burnin, n_thin = n_thin,
               frailty = dataL[[i]]$frailty)
  # get estimates
  jm_alphaF <- jm_fit$mcmc$alphaF
  jm_sigmaF <- jm_fit$mcmc$sigmaF
  jm_est <- list(alphaF = jm_alphaF, sigmaF = jm_sigmaF)
  # return
  list(jm_est = jm_est)
}
tic1 <- Sys.time()
set.seed(2021);seq_sample <- sample(seq_len(n_data), 20) 
cl_out <- parallel::makeCluster(ncl_out)
invisible(parallel::clusterEvalQ(cl_out, library(JMbayes2))) # load packages in each cl
parallel::clusterExport(cl_out, list("dataL", "scale", "n_iter", "n_burnin", "n_thin")) # load vars in each cl
res <- parallel::parLapply(cl_out, seq_sample, outer, ncl_in)
parallel::stopCluster(cl_out)
toc <- Sys.time()
saveRDS(res, file = paste0("Gen_Data/", ref, "/estimates5_", ref,".rds"))
RPushbullet::pbPost("note", title = paste0("RecData: fit complete"))
print(round(difftime(toc, tic1), 2))