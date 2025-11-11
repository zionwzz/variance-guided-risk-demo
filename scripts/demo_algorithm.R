library(ggplot2)
library(dplyr)
library(tidyr)

# ---- locate and source your code ----
repo_root <- getwd()
r_dir <- file.path(repo_root, "R")
stopifnot(dir.exists(r_dir),
          file.exists(file.path(r_dir, "simulator_scenarios.r")),
          file.exists(file.path(r_dir, "temporalVarGuid1.r")))
source(file.path(r_dir, "simulator_scenarios.r"))  # provides simulate_scenario_smallT()
source(file.path(r_dir, "temporalVarGuid1.r"))     # provides lmvt(), predict.lmvt()

# ---- simple arg parsing ----
args <- commandArgs(trailingOnly = TRUE)
scenario_id <- if (length(args) >= 1) as.integer(args[[1]]) else 2L
S            <- if (length(args) >= 2) as.integer(args[[2]]) else 5L
Tlen         <- if (length(args) >= 3) as.integer(args[[3]]) else 200L
seed         <- if (length(args) >= 4) as.integer(args[[4]]) else 42L
stopifnot(scenario_id %in% 1:4)

# ---- helpers ----
pstdt <- function(q, nu) { stopifnot(nu > 2); pt(q * sqrt(nu/(nu - 2)), df = nu) }
.get_innov_spec <- function(sim){
  cfg  <- try(attr(sim, "params")$cfg, silent = TRUE)
  dist <- try(cfg$innov$dist, silent = TRUE); if (inherits(dist, "try-error") || is.null(dist)) dist <- "norm"
  nu   <- try(cfg$innov$nu,   silent = TRUE); if (inherits(nu,   "try-error") || is.null(nu))   nu   <- NA_real_
  list(dist = if (identical(dist, "norm")) "norm" else "t", nu = nu)
}
.xcols <- function(sim) setdiff(names(sim), c("s","t","y","mu","sigma","pi_true","scenario"))

# ---- simulate ONE draw from the chosen scenario ----
set.seed(seed)
sim <- simulate_scenario_smallT(
  scen = scenario_id,
  S = S, T = Tlen,
  L = 20,             # typical regime run length
  d_noise = 100,      # nuisance covariates 
  noise_kind = "ar1",
  threshold_c = 1.0
)
spec <- .get_innov_spec(sim)

# ---- fit the variance-guided model ----
df <- sim[, c("s","t","y", .xcols(sim)), drop = FALSE]

fit <- lmvt(
  data = df,
  p = 1, q = 0, r = 1, s_ord = 1,     
  lambda_beta  = 0.01,                
  lambda_gamma = 0.01,                
  standardize_X = TRUE,
  use_x_in_variance = TRUE,           
  verbose = TRUE
)

pred <- predict.lmvt(
  object   = fit,
  newdata  = df,
  threshold = 0,
  innov_g  = identical(spec$dist, "norm"),
  innov_t  = identical(spec$dist, "t"),
  df_t     = if (is.na(spec$nu)) 6 else spec$nu
)

# ---- convert to exceedance risk at a single cutoff (global 90th percentile) ----
c_cut <- quantile(df$y, 0.90, na.rm = TRUE)
z <- (c_cut - pred$muhat) / pmax(pred$sigmahat, 1e-8)
p_hat <- if (identical(spec$dist, "norm")) 1 - pnorm(z) else 1 - pstdt(z, spec$nu)

demo <- df %>%
  mutate(muhat = pred$muhat,
         sigmahat = pred$sigmahat,
         p_hat = p_hat,
         exceed = as.integer(y > c_cut))

# ---- quick plots for a single subject (s=1) ----
out_dir <- file.path(repo_root, "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

s_show <- 1
d1 <- demo %>% filter(s == s_show)

p1 <- ggplot(d1, aes(t, y)) +
  geom_line(linewidth = 0.4) +
  geom_ribbon(aes(ymin = muhat - 2*sigmahat, ymax = muhat + 2*sigmahat), alpha = 0.15) +
  geom_line(aes(y = muhat), linetype = 2) +
  geom_hline(yintercept = c_cut, linetype = 3) +
  labs(
    title = sprintf("Scenario %d • Subject %d: y, μ̂ and ±2σ̂ band", scenario_id, s_show),
    x = "time", y = "value"
  ) +
  theme_minimal(base_size = 12)

p2 <- ggplot(d1, aes(t, p_hat)) +
  geom_line(linewidth = 0.4) +
  geom_point(aes(y = exceed), alpha = 0.35, size = 0.8) +
  labs(
    title = sprintf("Scenario %d • Subject %d: exceedance risk P(Y>c)", scenario_id, s_show),
    x = "time", y = "risk / indicator"
  ) +
  theme_minimal(base_size = 12)

f1 <- file.path(out_dir, sprintf("scen%d_subj%d_y_mu_band.png", scenario_id, s_show))
f2 <- file.path(out_dir, sprintf("scen%d_subj%d_risk_vs_hits.png", scenario_id, s_show))
ggsave(f1, p1, width = 8, height = 3.3, dpi = 150)
ggsave(f2, p2, width = 8, height = 3.3, dpi = 150)

# ---- quick console summary ----
cat(sprintf("\nDemo run complete.\nScenario = %d, S = %d, T = %d, seed = %d\n", scenario_id, S, Tlen, seed))
cat(sprintf("Cutoff c (90th %% of pooled y): %.3f\n", c_cut))
cat(sprintf("Subject %d: realized exceedance rate = %.3f; mean predicted risk = %.3f\n",
            s_show, mean(d1$exceed), mean(d1$p_hat)))
cat(sprintf("Figures saved:\n  %s\n  %s\n\n", f1, f2))
