ref <- "01"
desc <- "gap time"
################################################################################
## fit data
dataL <- readRDS(file = paste0("Gen_Data/dataL_", ref,".rds"))
n_data <- 2 #length(dataL)
n_cores <- 4 #max(parallel::detectCores() - 1L, 1L)
#
n_chains <- 3L
ncl_in <- min(n_chains, n_cores) # number of inner clusters (per outer cluster)
ncl_out <- max(floor(n_cores / ncl_in), 1) # number of outer clusters
outer <- function(i, ncl_in){
  tic2 <- Sys.time()
  lme_fit <- lme(y ~ ns(time, k =  c(1, 3), B = c(0, 7) ), 
                 random = list(id = pdDiag(form = ~ ns(time, k = c(1, 3), 
                                                       B = c(0, 7)))),
                 data = dataL[[i]]$long, 
                 control = lmeControl(opt = "optim", niterEM = 45))
  cox_fit <- coxph(Surv(tstart, tstop, status) ~ (group + age) * strata(strata),
                   data = dataL[[i]]$rec_strt)
  jm_fit <- jm(cox_fit, lme_fit, time_var = "time", 
               functional_forms = list("y" = ~ value(y) * strata),
               cores = ncl_in, recurrent = TRUE)
  timer <- difftime(Sys.time(), tic2, units = "mins")
  list(jm = jm_fit, timer = timer)
}
tic1 <- Sys.time()
cl_out <- parallel::makeCluster(ncl_out)
invisible(parallel::clusterEvalQ(cl_out, library(JMbayes2))) # load packages in each cl
parallel::clusterExport(cl_out, list("dataL")) # load vars in each cl
res <- parallel::parLapply(cl_out, seq_len(n_data), outer, ncl_in)
parallel::stopCluster(cl_out)
toc <- Sys.time()
diff_time <- sapply(res, "[[", "timer")
out <- lapply(res, "[[", "jm")
source("Gen_Data/sim_funs.R")
par_diag(n_data, n_chains, ncl_in, ncl_out, n_cores, diff_time)
print(round(difftime(toc, tic1), 2))
