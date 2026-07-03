# ============================================================================
# make_figures.R
# Regenerate the manuscript figures in the paper's ggplot2 style (six-model
# simulation boxplots; WESAD threshold-sensitivity and calibration). No titles.
# Run from the repository root:  Rscript manuscript/R/make_figures.R
# ============================================================================

ROOT <- normalizePath(".")
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

REP  <- file.path(ROOT, "manuscript", "data")
OUT  <- file.path(ROOT, "manuscript", "outputs")
SUPA <- file.path(ROOT, "manuscript", "outputs")
FIGM <- file.path(ROOT, "manuscript", "outputs", "figures")
FIGS <- file.path(ROOT, "manuscript", "outputs", "figures")
dir.create(FIGM, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)

## base-PDF device to avoid Type-3 fonts (as in the original manuscript code)
safe_pdf <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE, ...)

## ---- model labels & colours (A-D = original 4-hue; E,F added) --------------
model_labels <- c("A"="A: ARX–GARCH", "B"="B: ARXp–GARCH",
                  "C"="C: ARX–GARCHX", "D"="D: ARXp–GARCHXp",
                  "E"="E: ARXp–GARCHX", "F"="F: AR–GARCH")
model_cols   <- c("A"="#F8766D","B"="#7CAE00","C"="#00BFC4",
                  "D"="#C77CFF","E"="#E69F00","F"="#666666")

## ---- assemble per-replication long data (A-D published + E/F new) ----------
ef_raw <- read.csv(file.path(OUT, "sim_EF_raw.csv"), check.names = FALSE)
long_for <- function(metric) {              # metric = "RMSE" or "Bias"
  do.call(rbind, lapply(c(120, 240), function(Tt) {
    ad <- read.csv(file.path(REP, sprintf("rep_raw_s15t%d_by_rep_100reps.csv", Tt)),
                   check.names = FALSE)
    ef <- ef_raw[ef_raw$T == Tt, c("scenario_id","rep","E_RMSE","F_RMSE","E_Bias","F_Bias")]
    d  <- merge(ad, ef, by = c("scenario_id","rep"))
    cols <- paste0(c("A","B","C","D","E","F"), "_", metric)
    d %>% dplyr::select(scenario_id, rep, dplyr::all_of(cols)) %>%
      tidyr::pivot_longer(dplyr::all_of(cols), names_to = "mm", values_to = "value") %>%
      dplyr::mutate(model = factor(sub(paste0("_", metric), "", mm),
                                   levels = c("A","B","C","D","E","F")),
                    T = factor(paste0("T = ", Tt), levels = c("T = 120","T = 240")))
  }))
}

## ---- the shared plotting style (matches the original manuscript panels) -----
box_style <- function(df, ylab, ylim_use) {
  ggplot(df, aes(x = model, y = value, fill = model)) +
    geom_boxplot(outlier.alpha = 0.35, width = 0.7) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2.4, fill = "white") +
    facet_grid(T ~ scenario_id, labeller = labeller(scenario_id = function(s) paste0("S", s))) +
    scale_fill_manual(values = model_cols, labels = model_labels) +
    labs(x = "Model", y = ylab, fill = "Model") +           # NO title (added in LaTeX)
    coord_cartesian(ylim = ylim_use) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title      = element_text(size = 16),
      axis.text.x     = element_text(size = 11),
      axis.text.y     = element_text(size = 12),
      strip.text      = element_text(size = 13, face = "bold"),
      legend.title    = element_text(size = 13),
      legend.text     = element_text(size = 12),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

p_rmse <- box_style(long_for("RMSE"), "RMSE", c(0, 0.17))
p_bias <- box_style(long_for("Bias"), "Bias", c(-0.05, 0.05))
ggsave(file.path(FIGM, "Figure3_sixmodels_rmse.pdf"), p_rmse, width = 11, height = 6.4, units = "in", device = safe_pdf)
ggsave(file.path(FIGM, "Figure3_sixmodels_bias.pdf"), p_bias, width = 11, height = 6.4, units = "in", device = safe_pdf)
cat("Wrote six-model Figure 3 (bias, rmse) to", FIGM, "\n")

## ---- Supplement Fig S1: threshold-choice sensitivity (no title) ------------
ts <- read.csv(file.path(SUPA, "threshold_sensitivity.csv"), check.names = FALSE)
ord <- c("full-timeframe q0.70","burn-in first 20% q0.70","burn-in first 30% q0.70","baseline-window q0.70")
ts$threshold_rule <- factor(ts$threshold_rule, levels = ord)
ts$lab <- c("full timeframe\n(main analysis)","burn-in\nfirst 20%","burn-in\nfirst 30%","baseline\nwindows")[match(ts$threshold_rule, ord)]
ts$lab <- factor(ts$lab, levels = ts$lab[order(ts$threshold_rule)])
p_thr <- ggplot(ts, aes(x = lab, y = AUC_stress, fill = threshold_rule)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_hline(yintercept = 0.870, linetype = 2, colour = "grey40") +
  geom_text(aes(label = sprintf("%.3f", AUC_stress)), vjust = -0.5, size = 4) +
  annotate("text", x = 4.3, y = 0.872, label = "raw HR = 0.870", hjust = 1, size = 3.4, colour = "grey40") +
  scale_fill_manual(values = c("#2c7fb8","#7fcdbb","#7fcdbb","#fdae61")) +
  coord_cartesian(ylim = c(0.80, 0.96)) +
  labs(x = NULL, y = "AUC (stress vs non-stress)") +       # NO title
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 11))
ggsave(file.path(FIGS, "threshold_sensitivity.pdf"), p_thr, width = 6.6, height = 4.0, units = "in", device = safe_pdf)

## ---- Supplement Fig S2: Model-D reliability (no title) ---------------------
rel <- read.csv(file.path(SUPA, "modelD_calibration_reliability.csv"), check.names = FALSE)
p_cal <- ggplot(rel, aes(x = mean_pred, y = obs_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_line(colour = "#2c7fb8", linewidth = 0.7) +
  geom_point(colour = "#2c7fb8", size = 2.2) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Mean predicted risk (decile)", y = "Observed exceedance rate") +   # NO title
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FIGS, "modelD_calibration.pdf"), p_cal, width = 4.6, height = 4.6, units = "in", device = safe_pdf)

cat("Wrote Supplement figures S1 (threshold) and S2 (calibration) to", FIGS, "\n")
cat("Done. Recompile the manuscript and supplement to pick up the R-styled figures.\n")
