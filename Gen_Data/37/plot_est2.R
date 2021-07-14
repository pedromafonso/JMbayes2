setwd( "C:/Users/pedro/Documents/GitHub/JMbayes2-RE")
ref <- "37"
################################################################################
estimates <- readRDS(file = paste0("Gen_Data/", ref, "/estimates5_", ref,".rds"))
jm_est  <- lapply(estimates, "[[", "jm_est")

set.seed(2021);seq_sample <- sample(100, 20) 

# alphaF
alphaF_mcmc <- lapply(jm_est, "[[", "alphaF")
alphaF_true <- 0.25
par(mfrow = c(4, 5))
ylim <- range(unlist(alphaF_mcmc), alphaF_true)
xlim <- c(1, length(alphaF_mcmc[[1]][[1]]))
for(i in seq_along(alphaF_mcmc)) {
  plot(NA, type = "n", xlim = xlim, ylim = ylim, 
       xlab = "Iteration", ylab = "Estimate")
  lines(seq_along(alphaF_mcmc[[i]][[1]]), alphaF_mcmc[[i]][[1]], col = 2, lty = 2)
  lines(seq_along(alphaF_mcmc[[i]][[2]]), alphaF_mcmc[[i]][[2]], col = 3, lty = 2)
  lines(seq_along(alphaF_mcmc[[i]][[3]]), alphaF_mcmc[[i]][[3]], col = 4, lty = 2)
  abline(h = alphaF_true, lwd = 2)
  mtext(paste0("Dataset ", seq_sample[i]), side = 3, line = 0.5, cex = 0.7)
}
mtext("alphaF", side = 3, font = 2, line = -2, outer = TRUE)

# sigmaF
sigmaF_mcmc <- lapply(jm_est, "[[", "sigmaF")
sigmaF_true <- 0.25
par(mfrow = c(4, 5))
ylim <- range(unlist(sigmaF_mcmc), sigmaF_true)
xlim <- c(1, length(sigmaF_mcmc[[1]][[1]]))
for(i in seq_along(sigmaF_mcmc)) {
  plot(NA, type = "n", xlim = xlim, ylim = ylim, 
       xlab = "Iteration", ylab = "Estimate")
  lines(seq_along(sigmaF_mcmc[[i]][[1]]), sigmaF_mcmc[[i]][[1]], col = 2, lty = 2)
  lines(seq_along(sigmaF_mcmc[[i]][[2]]), sigmaF_mcmc[[i]][[2]], col = 3, lty = 2)
  lines(seq_along(sigmaF_mcmc[[i]][[3]]), sigmaF_mcmc[[i]][[3]], col = 4, lty = 2)
  abline(h = sigmaF_true, lwd = 2)
  mtext(paste0("Dataset ", seq_sample[i]), side = 3, line = 0.5, cex = 0.7)
}
mtext("sigmaF", side = 3, font = 2, line = -2, outer = TRUE)

