ref <- "01"
desc <- "gap time"
################################################################################
## generate data
n_data <- 4L
n_cores <- max(parallel::detectCores() - 1, 1L)
#
tic <- Sys.time()
fun <- function(i) {
  source("Gen_Data/sim_funs.R")
  seed <- 2021 + i
  data <- gen_RecData(seed = seed, n_scl = 1.5, 
                      alpha_r = 0.5, alpha_t = 0.5, scale = "gap")
}
cl <- parallel::makeCluster(n_cores)
dataL <- parallel::parLapply(cl, seq_len(n_data), fun)
parallel::stopCluster(cl)
toc <- Sys.time()
print(round(difftime(toc, tic), 2))
saveRDS(dataL, file = paste0("Gen_Data/dataL_", ref,".rds"))
