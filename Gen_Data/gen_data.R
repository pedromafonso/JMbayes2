setwd( "C:/Users/pedro/Documents/GitHub/JMbayes2-RE")
ref <- "11"
desc <- "valued frailty and alphaF"
n_data <- 100L
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
dir.create(paste0("Gen_Data/", ref))
saveRDS(dataL, file = paste0("Gen_Data/", ref, "/dataL_", ref,".rds"))
dur_min <- difftime(toc, tic, units = "min")
invisible(file.copy(from = "Gen_Data/fit_data.R", 
                    to = paste0("Gen_Data/", ref, "/gen_data.R"), 
                    overwrite = TRUE))
invisible(file.copy(from = "Gen_Data/sim_funs.R", 
                    to = paste0("Gen_Data/", ref, "/sim_funs.R"), 
                    overwrite = TRUE))
RPushbullet::pbPost("note", 
                    title = paste0("RecData: gen complete (", 
                                   round(dur_min, 2), " min)"))
print(round(difftime(toc, tic), 2))
