library(dplyr)
library(tidyr)
library(ggplot2)

repo_root <- getwd()
r_dir <- file.path(repo_root, "R")
stopifnot(dir.exists(r_dir),
          file.exists(file.path(r_dir, "simulator_scenarios.r")),
          file.exists(file.path(r_dir, "temporalVarGuid1.r")))
source(file.path(r_dir, "simulator_scenarios.r"))
source(file.path(r_dir, "temporalVarGuid1.r"))

# --------- helpers ---------
`%||%` <- function(x, y) if (is.null(x)) y else x
rmse     <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
pstdt    <- function(q, nu) { stopifnot(nu > 2); pt(q * sqrt(nu/(nu - 2)), df = nu) }
compute_thresholds <- function(y, probs = c(0.60,0.70,0.80,0.90,0.95), absolutes = NULL) {
  if (!is.null(absolutes)) {
    thr <- sort(unique(as.numeric(absolutes)))
    data.frame(thr_c = thr, thr_q = stats::ecdf(y)(thr), thr_type = "absolute")
  } else {
    thr <- as.numeric(stats::quantile(y, probs = probs, na.rm = TRUE, names = FALSE))
    data.frame(thr_c = thr, thr_q = probs, thr_type = "quantile")
  }
}
p_exceed <- function(mu, sigma, c, dist = "norm", nu = NA_real_) {
  z <- (c - mu) / pmax(sigma, 1e-8)
  if (identical(dist, "norm")) 1 - stats::pnorm(z) else 1 - pstdt(z, nu)
}
.get_innov_spec <- function(sim){
  cfg  <- try(attr(sim, "params")$cfg, silent = TRUE)
  dist <- try(cfg$innov$dist, silent = TRUE)
  nu   <- try(cfg$innov$nu,   silent = TRUE)
  if (inherits(dist, "try-error") || is.null(dist)) dist <- "norm"
  if (inherits(nu,   "try-error") || is.null(nu))   nu   <- NA_real_
  list(dist = if (identical(dist, "norm")) "norm" else "t", nu = nu)
}
.xcols <- function(sim) setdiff(names(sim), c("s","t","y","mu","sigma","pi_true","scenario"))
.ts    <- function() format(Sys.time(), "%H:%M:%S")

# --------- one model fit -> averaged RMSE & Bias over thresholds ---------
metrics_from_fit <- function(sim, use_x_in_variance,
                             lambda_beta, lambda_gamma,
                             p=1,q=0,r=1,s_ord=1,
                             standardize_X=TRUE,
                             use_sim_c = FALSE,
                             thr_probs = c(0.60,0.70,0.80,0.90,0.95),
                             verbose=FALSE) {
  spec  <- .get_innov_spec(sim)
  Xcols <- .xcols(sim)
  df <- sim[, c("s","t","y", Xcols), drop = FALSE]

  fit <- try(lmvt(
    data = df,
    p = p, q = q, r = r, s_ord = s_ord,
    lambda_beta  = lambda_beta,
    lambda_gamma = if (isTRUE(use_x_in_variance)) lambda_gamma else 0,
    standardize_X = standardize_X,
    use_x_in_variance = use_x_in_variance,
    verbose = verbose
  ), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(risk_rmse = Inf, risk_bias = Inf))

  pred <- try(predict.lmvt(
    object   = fit,
    newdata  = df,
    threshold = 0,
    innov_g  = identical(spec$dist, "norm"),
    innov_t  = identical(spec$dist, "t"),
    df_t     = if (is.na(spec$nu)) 6 else spec$nu
  ), silent = TRUE)
  if (inherits(pred, "try-error")) return(list(risk_rmse = Inf, risk_bias = Inf))

  thr_df <- if (isTRUE(use_sim_c)) {
    compute_thresholds(sim$y, absolutes = attr(sim, "params")$threshold %||% 1.0)
  } else {
    compute_thresholds(sim$y, probs = thr_probs)
  }

  rmse_vec <- bias_vec <- numeric(nrow(thr_df))
  for (k in seq_len(nrow(thr_df))) {
    cthr   <- thr_df$thr_c[k]
    p_true <- p_exceed(sim$mu,     sim$sigma,     c = cthr, dist = spec$dist, nu = spec$nu)
    p_hat  <- p_exceed(pred$muhat, pred$sigmahat, c = cthr, dist = spec$dist, nu = spec$nu)
    rmse_vec[k] <- rmse(p_hat, p_true)
    bias_vec[k] <- mean(p_hat - p_true, na.rm = TRUE)
  }

  list(risk_rmse = mean(rmse_vec), risk_bias = mean(bias_vec))
}

# --------- one replication over A/B/C/D ---------
evaluate_ABCD_fixed_once <- function(
  scen, seed,
  S=15, T=120, L=20, d_noise=100,
  lambda_B_beta = 0.01,
  lambda_D_beta = 0.01,
  lambda_D_gamma_by_scen = c(`1`=0.001,`2`=0.001,`3`=0.001,`4`=0.001),
  p=1,q=0,r=1,s_ord=1,
  standardize_X=TRUE,
  use_sim_c = FALSE,
  thr_probs=c(0.60,0.70,0.80,0.90,0.95),
  verbose=FALSE,
  noise_kind = c("ar1","iid")
){
  noise_kind <- match.arg(noise_kind)
  set.seed(seed)

  sim <- simulate_scenario_smallT(
    scen = scen, S = S, T = T, L = L, d_noise = d_noise,
    noise_kind = noise_kind, threshold_c = 1.0
  )
  attr(sim, "params")$cfg$scenario_id <- scen

  A <- metrics_from_fit(sim, FALSE, 0, 0, p,q,r,s_ord, standardize_X, use_sim_c, thr_probs, verbose)
  B <- metrics_from_fit(sim, FALSE, lambda_B_beta, 0, p,q,r,s_ord, standardize_X, use_sim_c, thr_probs, verbose)
  C <- metrics_from_fit(sim, TRUE,  0, 0, p,q,r,s_ord, standardize_X, use_sim_c, thr_probs, verbose)
  D <- metrics_from_fit(sim, TRUE,  lambda_D_beta,
                        lambda_D_gamma_by_scen[as.character(scen)],
                        p,q,r,s_ord, standardize_X, use_sim_c, thr_probs, verbose)

  data.frame(
    scenario_id = scen,
    A_RMSE = A$risk_rmse, B_RMSE = B$risk_rmse, C_RMSE = C$risk_rmse, D_RMSE = D$risk_rmse,
    A_Bias = A$risk_bias, B_Bias = B$risk_bias, C_Bias = C$risk_bias, D_Bias = D$risk_bias
  )
}

# --------- run replications (per T) ---------
run_for_T <- function(Tval,
                      n_reps = 30,
                      seed_base = 20251111,
                      S=15, L=20, d_noise=100,
                      thr_probs=c(0.60,0.70,0.80,0.90,0.95),
                      noise_kind = "ar1") {
  rows <- list()
  for (sc in 1:4) {
    cat(sprintf("[%s] T=%d scenario %d …\n", .ts(), Tval, sc)); flush.console()
    scen_rows <- vector("list", n_reps)
    for (rp in seq_len(n_reps)) {
      seed <- seed_base + sc*100000 + Tval*10 + rp
      scen_rows[[rp]] <- evaluate_ABCD_fixed_once(
        scen = sc, seed = seed,
        S=S, T=Tval, L=L, d_noise=d_noise,
        thr_probs = thr_probs,
        noise_kind = noise_kind
      )
    }
    rows[[as.character(sc)]] <- dplyr::bind_rows(scen_rows) %>% dplyr::mutate(rep = seq_len(dplyr::n()))
  }
  dplyr::bind_rows(rows) %>% dplyr::arrange(scenario_id, rep)
}

# ===== main run: T = 120 and 240 =====
out_dir <- file.path(repo_root, "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

res120 <- run_for_T(Tval = 120, n_reps = 30)
res240 <- run_for_T(Tval = 240, n_reps = 30)

# save raw
write.csv(res120, file.path(out_dir, "rep_raw_s15t120_by_rep_30reps.csv"), row.names = FALSE)
write.csv(res240, file.path(out_dir, "rep_raw_s15t240_by_rep_30reps.csv"), row.names = FALSE)

# ===== figures: boxplots (facets = scenarios; x = A/B/C/D) =====
plot_block <- function(df_raw, kind = c("RMSE","Bias"), Tval) {
  kind <- match.arg(kind)
  if (kind == "RMSE") {
    long <- df_raw %>%
      tidyr::pivot_longer(c(A_RMSE,B_RMSE,C_RMSE,D_RMSE), names_to="model_raw", values_to="val") %>%
      dplyr::mutate(model = factor(model_raw,
        levels=c("A_RMSE","B_RMSE","C_RMSE","D_RMSE"),
        labels=c("A","B","C","D")))
    ylab <- "Risk RMSE"
    title <- sprintf("Risk RMSE, S=15, T=%d, AR(1) covariate noise", Tval)
  } else {
    long <- df_raw %>%
      tidyr::pivot_longer(c(A_Bias,B_Bias,C_Bias,D_Bias), names_to="model_raw", values_to="val") %>%
      dplyr::mutate(model = factor(model_raw,
        levels=c("A_Bias","B_Bias","C_Bias","D_Bias"),
        labels=c("A","B","C","D")))
    ylab <- "Risk Bias"
    title <- sprintf("Risk Bias, S=15, T=%d, AR(1) covariate noise", Tval)
  }

  ggplot(long, aes(x=model, y=val, fill=model)) +
    geom_boxplot(outlier.alpha=.35, width=.7) +
    { if (kind == "Bias") geom_hline(yintercept=0, linetype=3) } +
    stat_summary(fun=median, geom="point", shape=23, size=2.2, fill="white") +
    labs(title=title, x="Model", y=ylab) +
    facet_wrap(~factor(scenario_id), nrow=1) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

g_rmse_120 <- plot_block(res120, "RMSE", 120)
g_bias_120 <- plot_block(res120, "Bias", 120)
g_rmse_240 <- plot_block(res240, "RMSE", 240)
g_bias_240 <- plot_block(res240, "Bias", 240)

ggsave(file.path(out_dir, "ABCD_s15_T120_ar1_RMSE.png"), g_rmse_120, width = 12, height = 3.2, dpi = 150)
ggsave(file.path(out_dir, "ABCD_s15_T120_ar1_Bias.png"), g_bias_120, width = 12, height = 3.2, dpi = 150)
ggsave(file.path(out_dir, "ABCD_s15_T240_ar1_RMSE.png"), g_rmse_240, width = 12, height = 3.2, dpi = 150)
ggsave(file.path(out_dir, "ABCD_s15_T240_ar1_Bias.png"), g_bias_240, width = 12, height = 3.2, dpi = 150)
