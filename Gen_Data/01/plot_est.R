ref <- "01"
desc <- "gap time"
n_data <- 200L
n <- 500L
scale <- "gap"
################################################################################
estimates <- readRDS(file = paste0("Gen_Data/estimates_", ref,".rds"))
dur_min   <- readRDS(file = paste0("Gen_Data/dur_min_", ref,".rds"))

lme_est <- lapply(estimates, "[[", "lme_est") 
jm_est  <- lapply(estimates, "[[", "jm_est")
lme_est <- do.call(rbind, lapply(lme_est, unlist))
colnames(lme_est) <- c("beta1", "beta2", "beta3", "beta4", 
                       "D[1,1]", "D[2,1]", "D[3,1]", "D[4,1]", "D[2,2]",
                       "D[2,3]", "D[2,4]", "D[3,3]", "D[3,4]", "D[4,4]", 
                       "sigma")
lme_est <- as.data.frame(lme_est)
lme_est$model <- "lme"
jm_est <- do.call(rbind, lapply(jm_est, unlist))
colnames(jm_est) <- c("beta1", "beta2", "beta3", "beta4",
                      "D[1,1]", "D[2,1]", "D[3,1]", "D[4,1]", "D[2,2]",
                      "D[2,3]", "D[2,4]", "D[3,3]", "D[3,4]", "D[4,4]", 
                      "sigma",
                      "alpha_R", "alpha_T", 
                      "gamma1_R", "gamma2_R", "gamma1_T", "gamma2_T")
jm_est[, "gamma1_T"] <- rowSums(jm_est[, c("gamma1_R", "gamma1_T")])
jm_est[, "gamma2_T"] <- rowSums(jm_est[, c("gamma2_R", "gamma2_T")])
jm_est[, "alpha_T"]  <- rowSums(jm_est[, c("alpha_R", "alpha_T")])

jm_est <- as.data.frame(jm_est)
jm_names <- names(jm_est)
jm_est$model <- "jm"

lme_est[, jm_names[!jm_names %in% names(lme_est)]] <- NA
prm_est <- rbind(jm_est, lme_est)

true_values <- c("beta1"    = 6.94,
                 "beta2"    = 1.30,
                 "beta3"    = 1.84,
                 "beta4"    = 1.82,
                 "sigma"    = 0.6,
                 "D[1,1]"   = 0.71, 
                 "D[2,1]"   = 0.33, 
                 "D[3,1]"   = 0.07, 
                 "D[4,1]"   = 1.26, 
                 "D[2,2]"   = 2.68,
                 "D[2,3]"   = 3.81, 
                 "D[2,4]"   = 4.35, 
                 "D[3,3]"   = 7.62, 
                 "D[3,4]"   = 5.4, 
                 "D[4,4]"   = 8,
                 "alpha_R"  = 0.5,
                 "alpha_T"  = 0.5, 
                 "gamma1_R" = 0.5,
                 "gamma2_R" = 0.05,
                 "gamma1_T" = 0.5,
                 "gamma2_T" = 0.05)

prms <- c("beta1", "beta2", "beta3", "beta4",
          "D[1,1]", "D[2,2]", "D[3,3]", "D[4,4]", 
          "sigma", "gamma1_R", "gamma2_R", "gamma1_T", "gamma2_T",
          "alpha_R", "alpha_T")

pdf(file = paste0("Gen_Data/plot_ref", ref,".pdf"),
    height = 10, width = 7.5)
{
  par(oma = c(2.5, 1, 1, 1))
  par(mfrow = c(4, 4), mar = c(2.1, 4.1, 4.1, 2.1))
  for(prm in prms) {
    ylim <- range(c(prm_est[, prm], true_values[prm]), na.rm = TRUE)
    boxplot(prm_est[, prm] ~ prm_est$model, ylim = ylim,
            main = prm, xlab = "", ylab = "", yaxt = "n")
    axis(2, las = 2, cex.axis = 0.75)
    abline(h = true_values[prm], col = 2)
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  n <- 500
  text <- c(paste0("Ref: ", ref, " (", format(Sys.Date(), format="%b %d"), ")"),
            paste0("Desc: ", desc),
            paste0("# datasets: ", n_data),
            paste0("# patients: ", n),
            paste0("Total dur (hr): ", round(dur_min/60, 2)),
            paste0("Dur/dt (min): ", round(dur_min/n_data, 2))
  )
  plot.new()
  legend("bottom", legend = text, bty = "n", ncol = 3)
}
dev.off()

