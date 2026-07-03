# ============================================================================
# 01_sim_add_EF.R
# Simulation study -- two further variants in the ARX-GARCH family:
#   E  ARXp-GARCHX : penalized mean, X-dependent (unpenalized) variance
#   F  AR-GARCH    : covariate-free baseline (no X in the mean or the variance)
# evaluated with the same generator, estimator, tuning, and seeds as Models A-D.
# Run from the repository root:  Rscript manuscript/R/01_sim_add_EF.R
# ============================================================================

ROOT <- normalizePath(".")
source(file.path(ROOT, "R", "simulator_scenarios.r"))  # simulate_scenario_smallT(), pstdt(), %||%
source(file.path(ROOT, "R", "temporalVarGuid.r"))     # lmvt(), predict.lmvt()
suppressMessages({library(glmnet); library(dplyr); library(tidyr)})
OUTD <- file.path(ROOT, "manuscript", "outputs")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

## ---- settings identical to the published A-D run --------------------------
SCEN        <- 1:4
T_GRID      <- c(120, 240)
N_REPS      <- 100               # FINAL RUN (validated at 5 reps: E≈D, F worst)
SEED_BASE   <- 20000422
S_SUBJ      <- 15
D_NOISE     <- 100
NOISE_KIND  <- "ar1"
THR_PROBS   <- c(0.60, 0.70, 0.80, 0.90, 0.95)
LAMBDA_D_BETA <- 0.01            # matches the published Model D mean penalty
USE_PARALLEL  <- TRUE

## ---- exact helpers copied verbatim from dissertation.qmd ------------------
rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
p_exceed <- function(mu, sigma, c, dist, nu) {
  z <- (c - mu) / pmax(sigma, 1e-8)
  if (identical(dist, "norm")) 1 - pnorm(z) else 1 - pstdt(z, nu = nu)
}
compute_thresholds <- function(y, probs = THR_PROBS, absolutes = NULL) {
  thr <- as.numeric(quantile(y, probs = probs, na.rm = TRUE, names = FALSE))
  data.frame(thr_c = thr, thr_q = probs, thr_type = "quantile")
}
xcols_from_sim <- function(df) setdiff(names(df), c("s","t","y","mu","sigma","pi_true","scenario"))
.xcols <- function(sim) xcols_from_sim(sim)
.get_innov_spec <- function(sim){
  cfg <- attr(sim, "params")$cfg
  list(dist = if (identical(cfg$innov$dist, "norm")) "norm" else "t",
       nu   = if (!is.null(cfg$innov$nu)) cfg$innov$nu else NA_real_)
}
metrics_from_fit <- function(sim, use_x_in_variance, lambda_beta, lambda_gamma,
                             p=1,q=0,r=1,s_ord=1, standardize_X=TRUE, thr_probs=THR_PROBS) {
  spec  <- .get_innov_spec(sim); Xcols <- .xcols(sim)
  df <- sim[, c("s","t","y", Xcols), drop = FALSE]
  fit <- try(lmvt(data = df, p=p,q=q,r=r,s_ord=s_ord,
                  lambda_beta = lambda_beta,
                  lambda_gamma = if (isTRUE(use_x_in_variance)) lambda_gamma else 0,
                  standardize_X = standardize_X, use_x_in_variance = use_x_in_variance), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(risk_rmse = Inf, risk_bias = Inf))
  pred <- try(predict.lmvt(object = fit, newdata = df, threshold = 0,
                  innov_g = identical(spec$dist, "norm"),
                  innov_t = identical(spec$dist, "t"),
                  df_t = if (is.na(spec$nu)) 6 else spec$nu), silent = TRUE)
  if (inherits(pred, "try-error")) return(list(risk_rmse = Inf, risk_bias = Inf))
  thr_df <- compute_thresholds(sim$y, probs = thr_probs)
  rr <- numeric(nrow(thr_df)); rb <- numeric(nrow(thr_df))
  for (k in seq_len(nrow(thr_df))) {
    cthr <- thr_df$thr_c[k]
    p_true <- p_exceed(sim$mu,     sim$sigma,     c = cthr, dist = spec$dist, nu = spec$nu)
    p_hat  <- p_exceed(pred$muhat, pred$sigmahat, c = cthr, dist = spec$dist, nu = spec$nu)
    rr[k] <- rmse(p_hat, p_true); rb[k] <- mean(p_hat - p_true, na.rm = TRUE)
  }
  list(risk_rmse = mean(rr), risk_bias = mean(rb))
}

## ---- E and F for one (scenario, T, rep) -----------------------------------
eval_EF_once <- function(scen, Tt, rep) {
  seed <- SEED_BASE + scen*100000 + rep     # identical scheme to the published run
  set.seed(seed)
  sim <- simulate_scenario_smallT(scen = scen, S = S_SUBJ, T = Tt, L = 20,
                                  d_noise = D_NOISE, noise_kind = NOISE_KIND, threshold_c = 1.0)
  attr(sim, "params")$cfg$scenario_id <- scen
  E  <- metrics_from_fit(sim, TRUE,  LAMBDA_D_BETA, 0)   # ARXp-GARCHX
  Fm <- metrics_from_fit(sim, FALSE, 1e3,           0)   # AR-GARCH
  data.frame(scenario_id = scen, T = Tt, rep = rep,
             E_RMSE = E$risk_rmse, F_RMSE = Fm$risk_rmse,
             E_Bias = E$risk_bias, F_Bias = Fm$risk_bias)
}

## ---- run (parallel if available) ------------------------------------------
tasks <- expand.grid(rep = seq_len(N_REPS), scen = SCEN, Tt = T_GRID)
run_k <- function(k) eval_EF_once(tasks$scen[k], tasks$Tt[k], tasks$rep[k])

have_fa <- USE_PARALLEL && requireNamespace("future.apply", quietly = TRUE) &&
           requireNamespace("future", quietly = TRUE)
if (have_fa) {
  future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))
  cat(sprintf("Running %d tasks in parallel on %d workers...\n",
              nrow(tasks), max(1, parallel::detectCores() - 1)))
  rows <- future.apply::future_lapply(seq_len(nrow(tasks)), run_k, future.seed = TRUE)
  future::plan(future::sequential)
} else {
  cat(sprintf("Running %d tasks serially...\n", nrow(tasks)))
  rows <- lapply(seq_len(nrow(tasks)), function(k){
    if (k %% 25 == 0) { cat("  task", k, "/", nrow(tasks), "\n"); flush.console() }
    run_k(k)
  })
}
ef <- dplyr::bind_rows(rows)
write.csv(ef, file.path(OUTD, "sim_EF_raw.csv"), row.names = FALSE)

## ---- summarise E/F (mean, se over reps) -----------------------------------
ef_summ <- ef %>%
  tidyr::pivot_longer(cols = matches("_(RMSE|Bias)$"),
                      names_to = c("model","metric"), names_sep = "_") %>%
  group_by(T, scenario_id, model, metric) %>%
  summarise(mean = mean(value), sd = sd(value), se = sd(value)/sqrt(n()), .groups = "drop")
write.csv(ef_summ, file.path(OUTD, "sim_EF_summary.csv"), row.names = FALSE)

cat("\n===== E/F RMSE mean by T x scenario =====\n")
print(as.data.frame(ef_summ %>% filter(metric == "RMSE") %>%
        select(T, scenario_id, model, mean) %>%
        tidyr::pivot_wider(names_from = model, values_from = mean)), digits = 4)
cat("\nWrote sim_EF_summary.csv and sim_EF_raw.csv to", OUTD, "\n")
cat("These E/F columns will be merged next to your published A-D in Table S1.\n")
