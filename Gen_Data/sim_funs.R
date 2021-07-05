# this code does not account (yet) for the presence of a frailty term
gen_RecData <- function (seed, n, n_scl = 1.5, alpha_r, alpha_t, scale = "gap") {
  
  set.seed(seed)
  n_target <- n # target number of observational units
  n <- n * n_scl
  n_i <- 15  # number of (planned) encounters per unit
  tmax <- 7 # maximum follow-up time (type I censoring)
  remove(n_scl)
  
  ##############################################################################
  # longitudinal outcome 1/2
  ## parameters true values
  betas <- c("Intercept" = 6.94, "Time1" = 1.30, "Time2" = 1.84, "Time3" = 1.82)
  sigma_y <- 0.6 # measurement error sd
  D <- matrix(0, 4, 4)
  D[lower.tri(D, TRUE)] <- c(0.71, 0.33, 0.07, 1.26, 2.68, 3.81, 4.35, 7.62, 5.4, 8)
  D <- D + t(D)
  diag(D) <- diag(D) * 0.5
  b <- MASS::mvrnorm(n, rep(0, nrow(D)), D)
  Bkn <- c(0, 7)
  kn <- c(1, 3)
  ## longitudinal data
  long <- data.frame(id   = rep(seq_len(n), each = n_i),
                     time = c(replicate(n, c(0, sort(runif(n_i - 1, 1, tmax))))))
  X <- model.matrix(~ 1 + splines::ns(time, knots = kn, Boundary.knots = Bkn), 
                    data = long)
  Z <- model.matrix(~ 1 + splines::ns(time, knots = kn, Boundary.knots = Bkn), 
                    data = long)
  eta_y <- as.vector(X %*% betas + rowSums(Z * b[long$id, ]))
  long$y <- rnorm(n * n_i, eta_y, sigma_y)
  remove(sigma_y, D, X, Z)
  
  ##############################################################################
  # terminal outcome
  ## parameters true values
  gammas_t <- c("(Intercept)" = -9, "Group" = 0.5, "Age" = 0.05) # phi = exp(Intercept)
  sigma_t <- 2
  alphaF <- 0.25 # association frailty
  sigmaF <- 0.25 # frailty SD
  frailty <- rnorm(n, mean = 0, sd = sigmaF)
  ## terminal data
  group <- rep(0:1, each = n/2)
  age <- runif(n, 30, 70)
  W_t <- cbind("(Intercept)" = 1, "Group" = group, "Age" = age)
  eta_t <- as.vector(W_t %*% gammas_t + alphaF * frailty) 
  invS_t <- function(t, u, i) {
    h <- function(s) { 
      NS <- splines::ns(s, knots = kn, Boundary.knots = Bkn)
      X <- cbind(1, NS)
      Z <- cbind(1, NS)
      eta_y <- as.vector(X %*% betas + rowSums(Z * b[rep(i, nrow(Z)), ]))
      exp(log(sigma_t) + (sigma_t - 1) * log(s) + eta_t[i] + eta_y * alpha_t) 
    }
    integrate(h, lower = 0, upper = t)$value + log(u)
  }
  u_t <- runif(n)
  ter_times <- numeric(n)
  for(i in seq_len(n)) {
    root <- try(uniroot(invS_t, interval = c(1e-05, 250), #?? update upper limit
                        u = u_t[i], i = i)$root, TRUE)  
    ter_times[i] <- if (!inherits(root, "try-error")) root else NA
  }
  surv_na <- !is.na(ter_times)
  if(sum(surv_na) < n_target) stop("Not enough patients. Increase 'n_scl'.")
  rmv_ids <- sample(which(surv_na), sum(surv_na) - n_target)
  surv_na[rmv_ids] <- FALSE # remove the excess of units
  surv <- data.frame(id    = seq_len(n)[surv_na],
                     time  = ter_times[surv_na],
                     group = group[surv_na],
                     age   = age[surv_na])
  cens_times <- tmax
  surv$Tstatus <- as.numeric(surv$time <= cens_times) # event indicator
  surv$time <- pmin(surv$time, cens_times) # add censoring time
  remove(gammas_t, sigma_t, group, W_t, eta_t, alpha_t, invS_t, u_t, i, root, n, 
         n_target, rmv_ids, ter_times, tmax, cens_times)
  
  ##############################################################################
  # longitudinal outcome 2/2
  long_na <- rep(surv_na, each = n_i)
  long <- long[long_na, , drop = FALSE] # drop patients removed in surv
  long_cens <- long$time <= rep(surv$time, each = n_i) 
  long <- long[long_cens, , drop = FALSE] # drop censored encounters
  remove(long_na, n_i, long_cens)
  
  ##############################################################################
  # recurring outcome
  ## parameters true values
  gammas_r <- c("(Intercept)" = -9+3, "Group" = 0.5, "Age" = 0.05) # phi = exp(Intercept)
  sigma_r <- 2
  ## recurring data
  W_r <- cbind("(Intercept)" = 1, "Group" = surv$group, "Age" = surv$age)
  eta_r <- as.vector(W_r %*% gammas_r + frailty)
  b <- b[surv_na, , drop = FALSE]
  if(scale == "gap") {
    invS_r <- function(t, u, i, tstart) {
      h <- function(s) { 
        # the difference is in the line below: s becomes s + tstart
        NS <- splines::ns(s + tstart, knots = kn, Boundary.knots = Bkn)
        X <- cbind(1, NS)
        Z <- cbind(1, NS)
        eta_y <- as.vector(X %*% betas + rowSums(Z * b[rep(i, nrow(Z)), ]))
        exp(log(sigma_r) + (sigma_r - 1) * log(s) + eta_r[i] + eta_y * alpha_r) 
      }
      integrate(h, lower = 0, upper = t)$value + log(u)
    }
  } else if(scale == "calendar") {
    invS_r <- function(t, u, i, tstart) {
      h <- function(s) { 
        # the difference is in the line below: s becomes s + tstart
        NS <- splines::ns(s + tstart, knots = kn, Boundary.knots = Bkn)
        X <- cbind(1, NS)
        Z <- cbind(1, NS)
        eta_y <- as.vector(X %*% betas + rowSums(Z * b[rep(i, nrow(Z)), ]))
        # the difference is in the line below: s becomes s + tstart
        exp(log(sigma_r) + (sigma_r - 1) * log(s + tstart) + eta_r[i] + eta_y * alpha_r) 
      }
      integrate(h, lower = 0, upper = t)$value + log(u)
    }
  }
  stop_times <- start_times <- id_times <- list()
  j <- 1
  for(i in seq_along(surv$id)) {
    tstart <- 0
    while(!is.na(tstart) & tstart < surv$time[i]) {
      u_r <- runif(1)
      root <- try(uniroot(invS_r, interval = c(1e-05, 250), #?? update upper limit
                          u = u_r, i = i, tstart = tstart)$root, TRUE)  
      tstop <- if(!inherits(root, "try-error")) root else NA
      start_times[[j]] <- tstart
      stop_times[[j]] <- tstart + tstop
      dur <- 0 #?? later add event duration (random or not-random)
      tstart <- tstart + tstop + dur
      id_times[[j]] <- surv$id[i]
      j <- j + 1
    }
  }
  rec <- data.frame(id     = unlist(id_times),                       
                    tstart = unlist(start_times),
                    tstop  = unlist(stop_times))
  
  long$id <- match(long$id, unique(long$id)) # rename IDs
  rec$id  <- match(rec$id, unique(rec$id))
  surv$id <- seq_along(surv$id)
  rec$group <- surv$group[rec$id]
  rec$age <- surv$age[rec$id]
  rec$Tstatus <- surv$Tstatus[rec$id]
  rec$Stime <- surv$time[rec$id]
  rec$Rstatus <- as.numeric(!is.na(rec$tstop) & rec$tstop < rec$Stime)  # event indicator
  rec$tstop <- pmin(rec$tstop, rec$Stime, na.rm = TRUE) # add cens time
  rec$gap <- rec$tstop - rec$tstart
  
  remove(gammas_r, sigma_r, W_r, eta_r, eta_y, alpha_r, surv_na, betas, b, 
         invS_r, stop_times, start_times, id_times, dur, j, i, tstart, u_r, 
         root, tstop, Bkn, kn)
  
  ##############################################################################
  # results to return
  #long$group <- surv$group[long$id]
  long$seed <- surv$seed <- rec$seed <- seed # save seed

  # rec & surv in strata format
  tail_rows <- cumsum(rle(rec$id)$length)
  new_rows <- sort(c(seq_along(rec$id), tail_rows))
  rec_strt <- rec[new_rows, ]
  rec_strt$strata <- 1
  rec_strt$status <- rec_strt$Rstatus
  tail_rows <- tail_rows + seq_along(tail_rows)
  rec_strt$strata[tail_rows] <- 2
  rec_strt$strata <- as.factor(rec_strt$strata)
  rec_strt$tstart[tail_rows] <- 0
  rec_strt$status[tail_rows] <- rec_strt$Tstatus[tail_rows]
  
  remove(seed, tail_rows, new_rows)
  rec_strt <- rec_strt[, c("id", "tstart", "tstop", "status", "strata", "Rstatus", "Tstatus", "Stime", "gap", "age", "group", "seed")]
  
  return(list(long = long, surv = surv, rec = rec, rec_strt = rec_strt, 
              frailty = frailty, alphaF = alphaF)) ##?? delete later
}

################################################################################

par_diag <- function(n_slc, n_chains, ncl_in, ncl_out, cores, diff_time) {
  # slice format
  nrow <- ceiling(n_chains / ncl_in)
  ncol <- ncl_in
  v <- seq_len(n_chains)
  length(v) <- nrow * ncol
  m <- matrix(v, byrow = TRUE, nrow = nrow, ncol = ncol)
  slc <- NULL
  for(r in seq_len(nrow)) {
    slc <- paste0(slc, m[r, !is.na(m[r, ])], collapse = "|")
    if(r < nrow) slc <- paste0(slc, "| \n ", collapse = "")
  }
  slc <- paste0("[", slc, "]")
  # diagram
  n_loops <- ceiling(n_slc / ncl_out)
  if(length(diff_time) != n_slc) {
    stop("The length of 'diff_time' doesn't match the expected size.")
  }
  length(diff_time) <- n_loops * ncl_out
  diff_time <- matrix(diff_time, byrow = TRUE, nrow = n_loops, ncol = ncl_out)
  diff_time <- apply(diff_time, 1, max, na.rm = TRUE)
  for(i in seq_len(n_loops)) {
    times <- min(n_slc - (i-1) * ncl_out, ncl_out)
    cat(rep(slc, times))
    free_cores <- cores - times * ncl_in
    if(free_cores > 0) {
      cat(paste0(" ", paste0(rep("_", free_cores), collapse = "|")))
    }
    cat(paste0(" ", round(diff_time[i], 2), "mins"))
    cat("\n\n")
  }
  
}

################################################################################