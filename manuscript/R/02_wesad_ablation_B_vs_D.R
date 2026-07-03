# ============================================================================
# 02_wesad_ablation_B_vs_D.R
# WESAD ablation: the full model (D: ARXp-GARCHXp, covariates in the variance)
# versus an otherwise identical model without covariates in the variance
# (B: ARXp-GARCH). Reports AUC and Brier for the exceedance event.
# Run from the repository root:  Rscript manuscript/R/02_wesad_ablation_B_vs_D.R
# ============================================================================

ROOT <- normalizePath(".")
EST  <- file.path(ROOT, "R", "temporalVarGuid.r")
DAT  <- file.path(ROOT, "manuscript", "data", "dat_fit.csv")       # 15 subjects x 77 features
LAB  <- file.path(ROOT, "manuscript", "data", "label_tbl.csv")     # s, t, label (2 = stress)
OUTD <- file.path(ROOT, "manuscript", "outputs")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

source(EST)
suppressMessages({library(dplyr); library(pROC)})

dat_fit   <- read.csv(DAT,  check.names = FALSE)
label_tbl <- read.csv(LAB,  check.names = FALSE)

## ---- fit the two models (same mean penalty; differ only in variance) ------
fitD <- lmvt(dat_fit, p = 1, q = 0, r = 1, s_ord = 1,
             lambda_beta = 0.01, lambda_gamma = 1e-4,
             use_x_in_variance = TRUE,  standardize_X = TRUE)   # X in variance
fitB <- lmvt(dat_fit, p = 1, q = 0, r = 1, s_ord = 1,
             lambda_beta = 0.01,
             use_x_in_variance = FALSE, standardize_X = TRUE)   # no X in variance

## ---- subject-specific 70th-percentile threshold (as in the paper) ---------
ids_chr <- as.character(fitD$ids)
y_by_s  <- split(dat_fit$y, dat_fit$s)
thr_s <- vapply(ids_chr, function(idc) as.numeric(quantile(y_by_s[[idc]], 0.70, na.rm = TRUE)),
                numeric(1)); names(thr_s) <- ids_chr

roc_best <- function(resp, score) {
  r1 <- pROC::roc(resp, score, levels = levels(resp), direction = "<", quiet = TRUE)
  r2 <- pROC::roc(resp, score, levels = levels(resp), direction = ">", quiet = TRUE)
  if (as.numeric(pROC::auc(r2)) >= as.numeric(pROC::auc(r1))) r2 else r1
}

evaluate <- function(fit, name) {
  # Get muhat/sigmahat with a scalar dummy threshold, then apply each subject's
  # own c_s(0.70) PER ROW (robust; avoids any per-subject-vector recycling).
  pr <- predict.lmvt(fit, newdata = dat_fit, threshold = 0, innov_g = TRUE, innov_t = FALSE)
  pr$cs   <- as.numeric(thr_s[as.character(pr$s)])
  pr$risk <- 1 - pnorm((pr$cs - pr$muhat) / pr$sigmahat)
  df <- pr %>% inner_join(label_tbl, by = c("s","t")) %>%
    inner_join(dplyr::select(dat_fit, s, t, y), by = c("s","t")) %>%
    filter(is.finite(risk)) %>%
    mutate(grp   = factor(ifelse(label == 2, "stress","non-stress"),
                          levels = c("non-stress","stress")),
           event = as.integer(y > cs))
  # (1) discrimination of the annotated stress state
  auc_stress <- as.numeric(pROC::auc(roc_best(df$grp, df$risk)))
  # (2) the ACTUAL exceedance event y > c_s(0.70): calibration + discrimination
  brier <- mean((df$risk - df$event)^2)
  auc_event <- as.numeric(pROC::auc(roc_best(
                 factor(ifelse(df$event==1,"stress","non-stress"),
                        levels=c("non-stress","stress")), df$risk)))
  data.frame(model = name, AUC_stress = round(auc_stress,3),
             AUC_event = round(auc_event,3), Brier_event = round(brier,4),
             event_rate = round(mean(df$event),3))
}

out <- rbind(evaluate(fitD, "D: ARXp-GARCHXp (X in variance)"),
             evaluate(fitB, "B: ARXp-GARCH  (no X in variance)"))
print(out)
write.csv(out, file.path(OUTD, "wesad_ablation_B_vs_D.csv"), row.names = FALSE)
cat("\nWrote", file.path(OUTD, "wesad_ablation_B_vs_D.csv"), "\n")
cat("Paste these AUC_stress / AUC_event / Brier_event values into the response and Supplementary Table.\n")
