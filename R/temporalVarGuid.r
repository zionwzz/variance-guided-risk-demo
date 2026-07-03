# ============================================================
# lmvt(): Penalized Panel ARX--GARCH--X with decoupled variance block
# ============================================================
# data: must contain columns s,t,y and any number of covariate columns (all non s,t,y are X)
# p,q   : mean lags (AR on y; X lags 0..q in the mean)
# r,s_ord: variance lags (ARCH/GARCH on eps^2 and sigma^2)
# lambda_beta: L1 on mean covariates (beta) via glmnet (if > 0); else weighted LS
# lambda_gamma: L1 on variance covariates (gamma) via proximal gradient (if > 0)
# standardize_X: standardize covariates (recommended for general use; FALSE for parameter recovery)
# use_x_in_variance: if FALSE, variance excludes X (gamma = 0, no B2 update)
# ------------------------------------------------------------

lmvt <- function(data,
                 p = 1, q = 0,            # mean lags
                 r = 1, s_ord = 1,        # variance lags
                 lambda_beta = 0, lambda_gamma = 0,
                 maxit = 200, tol = 1e-4,
                 standardize_X = TRUE,
                 use_x_in_variance = TRUE,
                 phi_cap = 0.995,
                 omega_min = 1e-6, omega_cap_mult = 10,
                 gamma_steps = 8,         # kept for compatibility
                 gamma_step = 1e-3,       # kept for compatibility
                 verbose = FALSE) {
  
  # ============================================================
  # If X is NOT used in variance, use original single-stage
  # ============================================================
  if (!use_x_in_variance) {
    return(lmvt_original_singlestage(data, p, q, r, s_ord, lambda_beta, lambda_gamma,
                                     maxit, tol, standardize_X, use_x_in_variance,
                                     phi_cap, omega_min, omega_cap_mult,
                                     gamma_steps, gamma_step, verbose))
  }
  
  # ============================================================
  # TWO-STAGE ESTIMATION (when use_x_in_variance = TRUE)
  # ============================================================
  if (verbose) cat("Using two-stage estimation for variance model with X...\n\n")
  
  # ============================================================
  # STAGE 1: Fit without X in variance
  # ============================================================
  if (verbose) cat("Stage 1: Estimating mean + ARCH/GARCH (no X in variance)...\n")
  
  fit_stage1 <- lmvt_original_singlestage(
    data = data, p = p, q = q, r = r, s_ord = s_ord,
    lambda_beta = lambda_beta,
    lambda_gamma = 0,
    use_x_in_variance = FALSE,  # KEY: no X in variance for stage 1
    standardize_X = standardize_X,
    maxit = maxit, tol = tol,
    phi_cap = phi_cap, omega_min = omega_min, omega_cap_mult = omega_cap_mult,
    gamma_steps = gamma_steps, gamma_step = gamma_step,
    verbose = verbose
  )
  
  if (verbose) {
    cat("\nStage 1 complete:\n")
    cat(sprintf("  theta = %.4f, beta[1] = %.4f\n", 
                if(length(fit_stage1$theta)) fit_stage1$theta[1] else 0,
                fit_stage1$beta[1]))
    cat(sprintf("  a = %.4f, b = %.4f, mean(omega) = %.4f\n\n",
                if(length(fit_stage1$a)) fit_stage1$a else 0,
                if(length(fit_stage1$b)) fit_stage1$b else 0,
                mean(fit_stage1$omega)))
  }
  
  # ============================================================
  # STAGE 2: Estimate gamma from residual variance
  # ============================================================
  if (verbose) cat("Stage 2: Estimating gamma from residual variance...\n")
  
  stopifnot(all(c("s","t","y") %in% names(data)))
  data <- data[order(data$s, data$t), ]
  ids <- fit_stage1$ids
  S <- length(ids)
  n <- nrow(data)
  
  # Reconstruct X
  xcols <- fit_stage1$xcols
  Xraw <- data.matrix(data[, xcols, drop = FALSE])
  if (fit_stage1$standardize_X) {
    X <- scale(Xraw, center = fit_stage1$x_center, scale = fit_stage1$x_scale)
  } else {
    X <- Xraw
  }
  Y <- as.numeric(data$y)
  
  # Lag helper
  lag_by_id <- function(v, L) {
    out <- rep(NA_real_, length(v))
    for (id in ids) {
      idx <- which(data$s == id)
      if (L <= 0) out[idx] <- v[idx] else if (length(idx) > L) {
        out[idx[(L+1):length(idx)]] <- v[idx[1:(length(idx)-L)]]
      }
    }
    out
  }
  
  # Create lags
  Ylags <- if (p > 0) do.call(cbind, lapply(1:p, function(L) lag_by_id(Y, L))) else NULL
  Xlags_list <- lapply(0:q, function(L) if (L == 0) X else apply(X, 2, lag_by_id, L = L))
  Xmean <- do.call(cbind, Xlags_list)
  Xvar <- Xmean
  
  # Get residuals from stage 1
  AR_part <- if (p > 0) as.numeric(Ylags %*% fit_stage1$theta) else 0
  X_part <- as.numeric(Xmean %*% fit_stage1$beta)
  mu_stage1 <- fit_stage1$alpha[match(data$s, ids)] + AR_part + X_part
  eps <- Y - mu_stage1
  eps[!is.finite(eps)] <- 0
  
  # Reconstruct sigma2 from stage 1
  sigma2_stage1 <- rep(NA_real_, n)
  
  make_var_lags <- function(eps, sigma2) {
    E2lags <- if (r > 0) do.call(cbind, lapply(1:r, function(L) lag_by_id(eps^2, L))) else NULL
    S2lags <- if (s_ord > 0) do.call(cbind, lapply(1:s_ord, function(L) lag_by_id(sigma2, L))) else NULL
    list(E2lags = E2lags, S2lags = S2lags)
  }
  
  for (id in ids) {
    idx <- which(data$s == id)
    for (ii in seq_along(idx)) {
      tt <- idx[ii]
      arch_part <- 0; garch_part <- 0
      if (r > 0) {
        for (L in 1:r) {
          tprev <- which(data$s[idx] == id & data$t[idx] == (data$t[tt] - L))
          if (length(tprev) == 1) {
            arch_part <- arch_part + fit_stage1$a[L] * (eps[idx[tprev]]^2)
          }
        }
      }
      if (s_ord > 0) {
        for (L in 1:s_ord) {
          tprev <- which(data$s[idx] == id & data$t[idx] == (data$t[tt] - L))
          if (length(tprev) == 1) {
            garch_part <- garch_part + fit_stage1$b[L] * sigma2_stage1[idx[tprev]]
          }
        }
      }
      sigma2_stage1[tt] <- max(fit_stage1$omega[match(id, ids)] + arch_part + garch_part, 1e-8)
    }
  }
  
  # Create mask for valid observations
  VL <- make_var_lags(eps, sigma2_stage1)
  mask_var <- rep(TRUE, n)
  if (r > 0) mask_var <- mask_var & (rowSums(is.na(VL$E2lags)) == 0)
  if (s_ord > 0) mask_var <- mask_var & (rowSums(is.na(VL$S2lags)) == 0)
  for (LL in 0:q) mask_var <- mask_var & (rowSums(is.na(Xlags_list[[LL+1]])) == 0)
  mask_var <- mask_var & is.finite(eps) & is.finite(sigma2_stage1)
  valid_idx <- which(mask_var)
  
  # Estimate gamma: regress (eps^2 - sigma2_stage1) on X
  residual_var <- (eps^2 - sigma2_stage1)[valid_idx]
  Xvar_valid <- Xvar[valid_idx, , drop = FALSE]
  
  # Weighted regression
  weights <- 1 / pmax(sigma2_stage1[valid_idx], 1e-6)
  weights <- weights / mean(weights)
  
  if (lambda_gamma > 0 && requireNamespace("glmnet", quietly = TRUE)) {
    if (verbose) cat(sprintf("  Using L1 penalty (lambda_gamma = %.4f)\n", lambda_gamma))
    Xw <- Xvar_valid * sqrt(weights)
    yw <- residual_var * sqrt(weights)
    fit_gamma <- glmnet::glmnet(Xw, yw, family = "gaussian", alpha = 1,
                                 lambda = lambda_gamma,
                                 standardize = FALSE, intercept = FALSE)
    gamma <- as.numeric(as.matrix(glmnet::coef.glmnet(fit_gamma)))[-1]
  } else {
    if (verbose) cat("  Using OLS\n")
    Xw <- Xvar_valid * sqrt(weights)
    yw <- residual_var * sqrt(weights)
    XtX <- crossprod(Xw)
    Xty <- crossprod(Xw, yw)
    gamma <- tryCatch(
      as.numeric(solve(XtX, Xty)),
      error = function(e) as.numeric(solve(XtX + 1e-6*diag(ncol(XtX)), Xty))
    )
  }
  
  gamma[!is.finite(gamma)] <- 0
  
  if (verbose) {
    cat("\nStage 2 complete:\n")
    cat(sprintf("  gamma[1:3] = %s\n",
                paste(round(gamma[1:min(3, length(gamma))], 4), collapse = ", ")))
    cat(sprintf("  Non-zero gammas: %d / %d\n\n", 
                sum(abs(gamma) > 1e-3), length(gamma)))
  }
  
  # Return model (same structure as original)
  structure(list(
    call = match.call(),
    orders = list(p = p, q = q, r = r, s = s_ord),
    ids = ids,
    alpha = fit_stage1$alpha,
    theta = fit_stage1$theta,
    beta = fit_stage1$beta,
    omega = fit_stage1$omega,
    a = fit_stage1$a,
    b = fit_stage1$b,
    gamma = gamma,
    xcols = xcols,
    standardize_X = fit_stage1$standardize_X,
    use_x_in_variance = TRUE,
    x_center = fit_stage1$x_center,
    x_scale = fit_stage1$x_scale,
    converged = fit_stage1$converged,
    iters = fit_stage1$iters
  ), class = "lmvt")
}


# ============================================================
# Helper: Original single-stage algorithm
# ============================================================
lmvt_original_singlestage <- function(data,
                 p = 1, q = 0,
                 r = 1, s_ord = 1,
                 lambda_beta = 0, lambda_gamma = 0,
                 maxit = 200, tol = 1e-4,
                 standardize_X = TRUE,
                 use_x_in_variance = TRUE,
                 phi_cap = 0.995,
                 omega_min = 1e-6, omega_cap_mult = 10,
                 gamma_steps = 8,
                 gamma_step = 1e-3,
                 verbose = FALSE) {
  
  stopifnot(all(c("s","t","y") %in% names(data)))
  data <- data[order(data$s, data$t), ]
  ids <- sort(unique(data$s)); S <- length(ids)
  n <- nrow(data)
  
  xcols <- setdiff(names(data), c("s","t","y"))
  if (length(xcols) == 0) stop("No covariates: provide columns beyond s,t,y.")
  Xraw <- data.matrix(data[, xcols, drop = FALSE])
  if (standardize_X) {
    x_center <- colMeans(Xraw, na.rm = TRUE)
    x_scale <- apply(Xraw, 2, sd, na.rm = TRUE)
    x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
    X <- scale(Xraw, center = x_center, scale = x_scale)
  } else {
    x_center <- rep(0, ncol(Xraw)); x_scale <- rep(1, ncol(Xraw)); X <- Xraw
  }
  P <- ncol(X)
  Y <- as.numeric(data$y)
  
  lag_by_id <- function(v, L) {
    out <- rep(NA_real_, length(v))
    for (id in ids) {
      idx <- which(data$s == id)
      if (L <= 0) out[idx] <- v[idx] else if (length(idx) > L) {
        out[idx[(L+1):length(idx)]] <- v[idx[1:(length(idx)-L)]]
      }
    }
    out
  }
  
  Ylags <- if (p > 0) do.call(cbind, lapply(1:p, function(L) lag_by_id(Y, L))) else NULL
  Xlags_list_mean <- lapply(0:q, function(L) if (L == 0) X else apply(X, 2, lag_by_id, L = L))
  Xmean <- do.call(cbind, Xlags_list_mean)
  
  if (use_x_in_variance) {
    Xlags_list_var <- Xlags_list_mean
    Xvar <- Xmean
    gamma <- rep(0, (q+1)*P)
  } else {
    Xlags_list_var <- list(matrix(0, nrow(data), 0))
    Xvar <- matrix(0, nrow(data), 0)
    gamma <- numeric(0)
    lambda_gamma <- 0
  }
  
  mask_mean <- rep(TRUE, n)
  if (p > 0) mask_mean <- mask_mean & (rowSums(is.na(Ylags)) == 0)
  for (LL in 0:q) mask_mean <- mask_mean & (rowSums(is.na(Xlags_list_mean[[LL+1]])) == 0)
  
  vl <- stats::var(Y, na.rm = TRUE); vl <- ifelse(is.finite(vl) && vl > 0, vl, 1)
  alpha <- tapply(Y, data$s, function(z) mean(z, na.rm = TRUE))
  alpha[is.na(alpha)] <- 0
  alpha <- as.numeric(alpha); names(alpha) <- ids
  omega <- rep(max(omega_min, vl/10), S); names(omega) <- ids
  theta <- if (p > 0) rep(0, p) else numeric(0)
  beta <- rep(0, (q+1)*P)
  a <- if (r > 0) rep(0.05, r) else numeric(0)
  b <- if (s_ord > 0) rep((phi_cap-0.1)/max(1, s_ord), s_ord) else numeric(0)
  
  mu <- alpha[match(data$s, ids)] +
    (if (p > 0) rowSums((Ylags*0), na.rm = TRUE) else 0) +
    as.numeric(Xmean %*% beta)
  eps <- Y - mu; eps[!is.finite(eps)] <- 0
  sigma2 <- rep(max(vl, 1e-2), n)
  
  make_var_lags <- function(eps, sigma2){
    E2lags <- if (r > 0) do.call(cbind, lapply(1:r, function(L) lag_by_id(eps^2, L))) else NULL
    S2lags <- if (s_ord > 0) do.call(cbind, lapply(1:s_ord, function(L) lag_by_id(sigma2, L))) else NULL
    list(E2lags = E2lags, S2lags = S2lags)
  }
  
  proj_ab <- function(a, b, cap = phi_cap) {
    a <- pmax(a, 0); b <- pmax(b, 0)
    s <- sum(a) + sum(b)
    if (s > cap && s > 0) { sc <- cap / s; a <- a*sc; b <- b*sc }
    list(a=a, b=b)
  }
  
  var_obj_abw <- function(par, eps, ids, s_idx, E2lags, S2lags, Xvar,
                          gamma, phi_cap, omega_min, omega_max, valid_idx) {
    BIG <- 1e50
    Sloc <- length(unique(ids))
    n_a <- if (is.null(E2lags)) 0 else ncol(E2lags)
    n_b <- if (is.null(S2lags)) 0 else ncol(S2lags)
    
    omega <- par[1:Sloc]
    a <- if (n_a > 0) par[(Sloc+1):(Sloc+n_a)] else numeric(0)
    b <- if (n_b > 0) par[(Sloc+n_a+1):(Sloc+n_a+n_b)] else numeric(0)
    
    if (any(!is.finite(omega))) return(BIG)
    omega <- pmin(pmax(omega, omega_min), omega_max)
    a[!is.finite(a)] <- 0; b[!is.finite(b)] <- 0
    if (any(a < 0) || any(b < 0)) return(BIG)
    ssum <- sum(a) + sum(b)
    if (ssum > phi_cap && ssum > 0) { sc <- phi_cap/ssum; a <- a*sc; b <- b*sc }
    
    delta <- omega[s_idx]
    if (n_a > 0) {
      if (anyNA(E2lags[valid_idx, , drop = FALSE])) return(BIG)
      delta <- delta + as.numeric(E2lags %*% a)
    }
    if (n_b > 0) {
      if (anyNA(S2lags[valid_idx, , drop = FALSE])) return(BIG)
      delta <- delta + as.numeric(S2lags %*% b)
    }
    u <- delta + if (length(gamma)) as.numeric(Xvar %*% gamma) else 0
    v <- valid_idx
    if (any(!is.finite(u[v]))) return(BIG)
    if (any(u[v] <= 0)) return(BIG)
    
    val <- sum(log(u[v]) + (eps[v]^2)/u[v])
    if (!is.finite(val)) return(BIG)
    val
  }
  
  obj_prev <- Inf
  for (k in 1:maxit) {
    
    w <- 1 / pmax(sigma2, 1e-10)
    AR_part <- if (p > 0) rowSums(Ylags %*% matrix(theta, ncol = 1)) else 0
    X_part <- as.numeric(Xmean %*% beta)
    numer <- tapply(w * (Y - AR_part - X_part), data$s, sum, na.rm = TRUE)
    denom <- tapply(w, data$s, sum, na.rm = TRUE)
    alpha <- as.numeric(numer / pmax(denom, 1e-12)); names(alpha) <- ids
    alpha[!is.finite(alpha)] <- 0
    
    if (p > 0) X_A <- cbind(Ylags, Xmean) else X_A <- Xmean
    y_A <- Y - alpha[match(data$s, ids)]
    keepA <- mask_mean & is.finite(y_A) & is.finite(w) &
      is.finite(rowSums(if (is.null(dim(X_A))) cbind(X_A) else X_A))
    X_Ak <- X_A[keepA, , drop = FALSE]; y_Ak <- y_A[keepA]; wk <- w[keepA]
    Xw <- X_Ak * sqrt(wk); yw <- y_Ak * sqrt(wk)
    
    if (lambda_beta > 0) {
      if (!requireNamespace("glmnet", quietly = TRUE))
        stop("Install 'glmnet' for penalized mean block.")
      pen_vec <- c(rep(0, p), rep(1, ncol(Xmean)))
      gfit <- glmnet::glmnet(Xw, yw, family="gaussian", alpha=1, lambda=lambda_beta,
                             penalty.factor=pen_vec, standardize=FALSE, intercept=FALSE)
      coef_all <- as.numeric(as.matrix(glmnet::coef.glmnet(gfit)))[-1]
    } else {
      XtX <- crossprod(Xw); Xty <- crossprod(Xw, yw)
      coef_all <- tryCatch(as.numeric(solve(XtX, Xty)),
                           error=function(e) as.numeric(solve(XtX + 1e-8*diag(ncol(XtX)), Xty)))
    }
    if (p > 0) { theta <- coef_all[1:p]; beta <- coef_all[(p+1):length(coef_all)] } else beta <- coef_all
    
    AR_part <- if (p > 0) rowSums(Ylags %*% matrix(theta, ncol = 1)) else 0
    X_part <- as.numeric(Xmean %*% beta)
    mu <- alpha[match(data$s, ids)] + AR_part + X_part
    eps <- Y - mu; eps[!is.finite(eps)] <- 0
    
    VL <- make_var_lags(eps, sigma2); E2lags <- VL$E2lags; S2lags <- VL$S2lags
    mask_var <- rep(TRUE, n)
    if (r > 0) mask_var <- mask_var & (rowSums(is.na(E2lags)) == 0)
    if (s_ord > 0) mask_var <- mask_var & (rowSums(is.na(S2lags)) == 0)
    if (use_x_in_variance) {
      for (LL in 0:q) mask_var <- mask_var & (rowSums(is.na(Xlags_list_var[[LL+1]])) == 0)
    }
    valid_idx <- which(mask_var & is.finite(eps))
    if (!length(valid_idx)) stop("No valid rows for variance update.")
    
    s_idx <- match(data$s, ids)
    omega_cap <- max(omega_cap_mult * stats::median(eps^2, na.rm = TRUE), omega_min*10)
    par0 <- c(omega, if (r>0) a, if (s_ord>0) b)
    lower <- c(rep(omega_min, S), if (r>0) rep(0, r), if (s_ord>0) rep(0, s_ord))
    upper <- c(rep(omega_cap, S), if (r>0) rep(phi_cap, r), if (s_ord>0) rep(phi_cap, s_ord))
    
    opt <- optim(par0, var_obj_abw, method = "L-BFGS-B",
                 lower = lower, upper = upper,
                 control = list(maxit = 200),
                 eps = eps, ids = data$s, s_idx = s_idx,
                 E2lags = E2lags, S2lags = S2lags, Xvar = Xvar,
                 gamma = gamma, phi_cap = phi_cap,
                 omega_min = omega_min, omega_max = omega_cap,
                 valid_idx = valid_idx)
    
    par_star <- opt$par
    omega <- par_star[1:S]
    if (r>0) a <- par_star[(S+1):(S+r)]
    if (s_ord>0) { off <- S + if (r>0) r else 0; b <- par_star[(off+1):(off+s_ord)] }
    abp <- proj_ab(a, b, phi_cap); a <- abp$a; b <- abp$b
    
    delta <- omega[s_idx] +
      if (r>0) as.numeric(E2lags %*% a) else 0 +
      if (s_ord>0) as.numeric(S2lags %*% b) else 0
    
    if (use_x_in_variance && length(gamma)) {
      for (gi in 1:max(1, gamma_steps)) {
        u <- delta + as.numeric(Xvar %*% gamma)
        v <- valid_idx
        uv <- u[v]
        step <- gamma_step
        for (bt in 0:8) {
          if (any(!is.finite(uv)) || any(uv <= 0)) { 
            step <- step * 0.5
            uv <- delta[v] + as.numeric(Xvar[v,] %*% gamma)
            next 
          }
          g_u <- (1/uv - (eps[v]^2)/(uv^2))
          grad_gam <- as.numeric(crossprod(Xvar[v, , drop=FALSE], g_u))
          gamma_try <- gamma - step * grad_gam
          if (lambda_gamma > 0) gamma_try <- sign(gamma_try) * pmax(abs(gamma_try) - step*lambda_gamma, 0)
          u_try <- delta[v] + as.numeric(Xvar[v, ] %*% gamma_try)
          if (any(!is.finite(u_try)) || any(u_try <= 0)) { step <- step * 0.5; next }
          obj_now <- sum(log(uv) + (eps[v]^2)/uv)
          obj_try <- sum(log(u_try) + (eps[v]^2)/u_try)
          if (is.finite(obj_try) && obj_try <= obj_now) { gamma <- gamma_try; break } else { step <- step * 0.5 }
        }
      }
    } else {
      gamma <- numeric(0)
    }
    
    sigma2_new <- rep(NA_real_, n)
    for (id in ids) {
      idx <- which(data$s == id)
      for (ii in seq_along(idx)) {
        tt <- idx[ii]
        arch_part <- 0; garch_part <- 0
        if (r > 0) {
          for (L in 1:r) {
            tprev <- which(idx == (data$t[tt] - L))
            if (length(tprev) == 1) arch_part <- arch_part + a[L] * (eps[idx[tprev]]^2)
          }
        }
        if (s_ord > 0) {
          for (L in 1:s_ord) {
            tprev <- which(idx == (data$t[tt] - L))
            if (length(tprev) == 1) garch_part <- garch_part + b[L] * sigma2_new[idx[tprev]]
          }
        }
        xrow_var <- if (ncol(Xvar)) Xvar[tt, ] else numeric(0)
        vlin <- omega[match(id, ids)] + arch_part + garch_part +
          if (length(gamma)) sum(gamma * xrow_var) else 0
        sigma2_new[tt] <- max(vlin, 1e-8)
      }
    }
    sigma2 <- sigma2_new
    
    u_now <- delta[valid_idx] + if (length(gamma)) as.numeric(Xvar[valid_idx, ] %*% gamma) else 0
    obj_now <- sum(log(u_now) + (eps[valid_idx]^2)/u_now)
    if (verbose) cat(sprintf("iter %d: obj=%s  sum(a)+sum(b)=%.4f  mean(omega)=%.4f\n",
                             k, format(obj_now, digits=6), sum(a)+sum(b), mean(omega)))
    if (!is.finite(obj_now)) break
    if (abs(obj_prev - obj_now) < tol) break
    obj_prev <- obj_now
  }
  
  structure(list(
    call = match.call(),
    orders = list(p=p, q=q, r=r, s=s_ord),
    ids = ids,
    alpha = alpha, theta = theta, beta = beta,
    omega = omega, a = a, b = b, gamma = gamma,
    xcols = xcols,
    standardize_X = standardize_X,
    use_x_in_variance = use_x_in_variance,
    x_center = x_center, x_scale = x_scale,
    converged = (k < maxit), iters = k
  ), class = "lmvt")
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
  AR_part <- if (p > 0) rowSums(Ylags %*% matrix(object$theta, ncol=1), na.rm=TRUE) else 0
  X_part  <- as.numeric(Xmean %*% object$beta)
  muhat   <- object$alpha[match(newdata$s, ids)] + AR_part + X_part
  eps_hat <- Y - muhat; eps_hat[!is.finite(eps_hat)] <- 0
  
  # forward variance recursion
  sigma2 <- rep(NA_real_, nrow(newdata))
  for (id in ids) {
    idx <- which(newdata$s == id)
    for (ii in seq_along(idx)) {
      tt <- idx[ii]
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
      xrow_var <- if (ncol(Xvar)) Xvar[tt, ] else numeric(0)
      vlin <- object$omega[match(id, ids)] + arch_part + garch_part +
        if (length(object$gamma)) sum(object$gamma * xrow_var) else 0
      sigma2[tt] <- max(vlin, 1e-8)
    }
  }
  
  sigmahat <- sqrt(sigma2); yhat <- muhat
  z <- (threshold - muhat)/sigmahat
  out <- data.frame(s = newdata$s, t = newdata$t,
                    yhat = yhat, muhat = muhat, sigmahat = sigmahat)
  if (innov_g) out$risk_g <- 1 - pnorm(z)
  if (innov_t) out$risk_t <- 1 - pt(z, df = df_t)
  out
}
