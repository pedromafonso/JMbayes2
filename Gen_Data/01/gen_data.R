ref <- "01"
desc <- "gap time"
n_data <- 200L
n <- 500L
scale <- "gap"
################################################################################
## generate data
n_cores <- max(parallel::detectCores() - 1, 1L)
tic <- Sys.time()
fun <- function(i) {
  source("Gen_Data/sim_funs.R")
  seed <- 2021 + i
  data <- gen_RecData(seed = seed, n = n, n_scl = 1.5, 
                      alpha_r = 0.5, alpha_t = 0.5, scale = scale)
}
cl <- parallel::makeCluster(n_cores)
parallel::clusterExport(cl, list("n", "scale")) # load vars in each cl
dataL <- parallel::parLapply(cl, seq_len(n_data), fun)
parallel::stopCluster(cl)
toc <- Sys.time()
saveRDS(dataL, file = paste0("Gen_Data/dataL_", ref,".rds"))
dur_min <- difftime(toc, tic, units = "min")
RPushbullet::pbPost("note", 
                    title = paste0("RecData: gen complete (", 
                                   round(dur_min, 2), " min)"))
print(round(difftime(toc, tic), 2))
