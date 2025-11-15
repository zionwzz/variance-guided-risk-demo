lmvt <- function(data,
                 p = 1, q = 0,            # mean lags
                 r = 1, s_ord = 1,        # variance lags
                 lambda_beta = 0, lambda_gamma = 0,
                 maxit = 200, tol = 1e-4,
                 standardize_X = TRUE,
                 use_x_in_variance = TRUE,  # NEW
                 phi_cap = 0.995,           # cap on sum(a)+sum(b)
                 omega_min = 1e-6, omega_cap_mult = 10,
                 gamma_steps = 8,           # inner prox steps per outer iter
                 gamma_step = 1e-3,         # prox step size for gamma
                 verbose = FALSE) {
  
  # ---- basic checks / setup (safe to publish) ----
  stopifnot(all(c("s", "t", "y") %in% names(data)))
  data <- data[order(data$s, data$t), ]
  ids  <- sort(unique(data$s))
  S    <- length(ids)
  n    <- nrow(data)
  
  # autodetect and (optionally) standardize X
  xcols <- setdiff(names(data), c("s", "t", "y"))
  if (length(xcols) == 0L) {
    stop("No covariates: provide columns beyond s,t,y.")
  }
  Xraw <- data.matrix(data[, xcols, drop = FALSE])
  
  if (standardize_X) {
    x_center <- colMeans(Xraw, na.rm = TRUE)
    x_scale  <- apply(Xraw,  2, sd, na.rm = TRUE)
    x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
    X <- scale(Xraw, center = x_center, scale = x_scale)
  } else {
    x_center <- rep(0, ncol(Xraw))
    x_scale  <- rep(1, ncol(Xraw))
    X <- Xraw
  }
  
  P <- ncol(X)
  Y <- as.numeric(data$y)
  
  # simple helper (non-core)
  lag_by_id <- function(v, L) {
    out <- rep(NA_real_, length(v))
    for (id in ids) {
      idx <- which(data$s == id)
      if (L <= 0) {
        out[idx] <- v[idx]
      } else if (length(idx) > L) {
        out[idx[(L+1):length(idx)]] <- v[idx[1:(length(idx)-L)]]
      }
    }
    out
  }

# ============================================================
# predict.lmvt(): conditional mean, variance, exceedance risks
# ============================================================

predict.lmvt <- function(object, newdata, threshold,
                         innov_g = TRUE, innov_t = TRUE, df_t = 6) {
  stopifnot(inherits(object, "lmvt"))
  stopifnot(all(c("s","t","y") %in% names(newdata)))
  newdata <- newdata[order(newdata$s, newdata$t), ]
  ids <- object$ids
  
  # guard: unseen subjects
  unseen <- setdiff(unique(newdata$s), ids)
  if (length(unseen)) {
    stop(sprintf("predict.lmvt: newdata contains unseen subject IDs: %s",
                 paste(unseen, collapse = ", ")))
  }
  
  # ---- resolve thresholds to per-row vector ----
  resolve_threshold <- function(threshold, new_s, train_ids) {
    n <- length(new_s)
    # 1) single fixed threshold
    if (length(threshold) == 1L) {
      thr <- rep(as.numeric(threshold), n)
      if (!all(is.finite(thr))) stop("Non-finite threshold.")
      return(thr)
    }
    # 2) per-row vector (length == nrow(newdata))
    if (length(threshold) == n) {
      thr <- as.numeric(threshold)
      if (any(!is.finite(thr))) stop("Non-finite values in threshold vector.")
      return(thr)
    }
    # 3) per-subject vector (preferred). Use names if provided; otherwise infer order.
    if (!is.null(names(threshold))) {
      # treat names as character keys (robust to numeric IDs stored as chars)
      thr_map <- setNames(as.numeric(threshold), as.character(names(threshold)))
      thr <- thr_map[as.character(new_s)]
      if (any(is.na(thr))) {
        missing_ids <- unique(new_s[is.na(thr)])
        stop(sprintf("Thresholds missing for subjects: %s",
                     paste(missing_ids, collapse = ", ")))
      }
      if (any(!is.finite(thr))) stop("Non-finite values in named threshold vector.")
      return(thr)
    } else {
      # no names: try matching by subject set size
      new_ids_unique <- unique(new_s)
      if (length(threshold) == length(new_ids_unique)) {
        thr_map <- setNames(as.numeric(threshold), as.character(new_ids_unique))
        thr <- thr_map[as.character(new_s)]
        if (any(!is.finite(thr))) stop("Non-finite values in per-subject threshold vector.")
        return(thr)
      }
      if (length(threshold) == length(train_ids)) {
        thr_map <- setNames(as.numeric(threshold), as.character(train_ids))
        thr <- thr_map[as.character(new_s)]
        if (any(is.na(thr))) {
          missing_ids <- unique(new_s[is.na(thr)])
          stop(sprintf("Thresholds (aligned to training IDs) missing for subjects: %s",
                       paste(missing_ids, collapse = ", ")))
        }
        if (any(!is.finite(thr))) stop("Non-finite values in per-subject threshold vector.")
        return(thr)
      }
      stop("`threshold` must be: length 1; length nrow(newdata); or a per-subject vector whose length equals the number of subjects (in newdata or in training). Provide names(threshold) = subject IDs to be explicit.")
    }
  }
  
  thr_per_row <- resolve_threshold(threshold, newdata$s, ids)
  
  # align covariates exactly like training
  Xraw <- matrix(0, nrow = nrow(newdata), ncol = length(object$xcols))
  colnames(Xraw) <- object$xcols
  present <- intersect(object$xcols, names(newdata))
  if (length(present) > 0)
    Xraw[, match(present, object$xcols)] <- data.matrix(newdata[, present, drop = FALSE])
  X <- if (isTRUE(object$standardize_X)) {
    scale(Xraw, center = object$x_center, scale = object$x_scale)
  } else Xraw
  
  Y <- as.numeric(newdata$y)
  
  # helper
  lag_by_id <- function(v, L) {
    out <- rep(NA_real_, length(v))
    for (id in ids) {
      idx <- which(newdata$s == id)
      if (L <= 0) out[idx] <- v[idx] else if (length(idx) > L) {
        out[idx[(L+1):length(idx)]] <- v[idx[1:(length(idx)-L)]]
      }
    }
    out
  }
  
  p <- object$orders$p; q <- object$orders$q
  r <- object$orders$r; s_ord <- object$orders$s
  
  # mean-side lags
  Ylags <- if (p > 0) do.call(cbind, lapply(1:p, function(L) lag_by_id(Y, L))) else NULL
  Xlags_list_mean <- lapply(0:q, function(L) if (L == 0) X else apply(X, 2, lag_by_id, L = L))
  Xmean <- do.call(cbind, Xlags_list_mean)
  
  # variance-side design follows training choice
  if (isTRUE(object$use_x_in_variance)) {
    Xlags_list_var <- Xlags_list_mean
    Xvar <- Xmean
  } else {
    Xlags_list_var <- list(matrix(0, nrow(newdata), 0))
    Xvar <- matrix(0, nrow(newdata), 0)
  }
  
  # conditional mean & residuals
  AR_part <- if (p > 0) rowSums(Ylags %*% matrix(object$theta, ncol = 1), na.rm = TRUE) else 0
  X_part  <- as.numeric(Xmean %*% object$beta)
  muhat   <- object$alpha[match(newdata$s, ids)] + AR_part + X_part
  eps_hat <- Y - muhat; eps_hat[!is.finite(eps_hat)] <- 0
  
  # forward variance recursion with feasibility repair
  u_floor <- 1e-8
  sigma2 <- rep(NA_real_, nrow(newdata))
  for (id in ids) {
    idx <- which(newdata$s == id)
    for (ii in seq_along(idx)) {
      tt <- idx[ii]
      
      # ARCH/GARCH parts using available lags within subject
      arch_part <- 0; garch_part <- 0
      if (r > 0) {
        for (L in 1:r) {
          tprev <- which(idx == (newdata$t[tt] - L))
          if (length(tprev) == 1) arch_part <- arch_part + object$a[L] * (eps_hat[idx[tprev]]^2)
        }
      }
      if (s_ord > 0) {
        for (L in 1:s_ord) {
          tprev <- which(idx == (newdata$t[tt] - L))
          if (length(tprev) == 1) garch_part <- garch_part + object$b[L] * sigma2[idx[tprev]]
        }
      }
      
      # variance linear predictor
      u0 <- object$omega[match(id, ids)] + arch_part + garch_part
      xrow_var <- if (ncol(Xvar)) Xvar[tt, ] else numeric(0)
      xg <- if (length(object$gamma)) sum(object$gamma * xrow_var) else 0
      vlin <- u0 + xg
      
      # Feasibility repair (prediction-time analogue of training B2 guards)
      if (!is.finite(vlin) || vlin <= u_floor) {
        if (is.finite(xg) && xg < 0) {
          s <- (u_floor - u0) / xg  # xg < 0
          s <- min(1, max(0, s))
          vlin_try <- u0 + s * xg
          if (is.finite(vlin_try) && vlin_try > u_floor) {
            vlin <- vlin_try
          } else {
            vlin <- max(u0, u_floor)
          }
        } else {
          vlin <- max(u0, u_floor)
        }
      }
      
      sigma2[tt] <- vlin
    }
  }
  
  sigmahat <- sqrt(pmax(sigma2, u_floor))
  yhat <- muhat
  z <- (thr_per_row - muhat) / sigmahat
  
  out <- data.frame(s = newdata$s, t = newdata$t,
                    yhat = yhat, muhat = muhat, sigmahat = sigmahat,
                    threshold_used = thr_per_row)
  if (innov_g) out$risk_g <- 1 - pnorm(z)
  if (innov_t) out$risk_t <- 1 - pt(z, df = df_t)
  out
}
