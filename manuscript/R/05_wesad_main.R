# ============================================================================
# 05_wesad_main.R
# Main WESAD analysis (headline result): the variance-aware
# exceedance risk score (Model D: ARXp-GARCHXp) versus the raw 60-s HR_mean
# baseline, evaluated on the annotated stress state using each subject's own
# 70th-percentile threshold c_s(0.70).
# Run from the repository root:  Rscript manuscript/R/05_wesad_main.R
# ============================================================================

ROOT <- normalizePath(".")
EST  <- file.path(ROOT, "R", "temporalVarGuid.r")
DAT  <- file.path(ROOT, "manuscript", "data", "dat_fit.csv")      # 15 subjects x 77 features
LAB  <- file.path(ROOT, "manuscript", "data", "label_tbl.csv")    # s, t, label (2 = stress)
OUTD <- file.path(ROOT, "manuscript", "outputs")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

source(EST)
suppressMessages({library(dplyr); library(pROC)})

dat_fit   <- read.csv(DAT, check.names = FALSE)
label_tbl <- read.csv(LAB, check.names = FALSE)

## ---- fit the main variance-aware risk model (Model D) ---------------------
fitD <- lmvt(dat_fit, p = 1, q = 0, r = 1, s_ord = 1,
             lambda_beta = 0.01, lambda_gamma = 1e-4,
             use_x_in_variance = TRUE, standardize_X = TRUE)

## ---- subject-specific 70th-percentile threshold c_s(0.70) -----------------
ids_chr <- as.character(fitD$ids)
y_by_s  <- split(dat_fit$y, dat_fit$s)
thr_s <- vapply(ids_chr, function(idc) as.numeric(quantile(y_by_s[[idc]], 0.70, na.rm = TRUE)),
                numeric(1)); names(thr_s) <- ids_chr

## auto-orient ROC (score may run either way) --------------------------------
roc_best <- function(resp, score) {
  r1 <- pROC::roc(resp, score, levels = levels(resp), direction = "<", quiet = TRUE)
  r2 <- pROC::roc(resp, score, levels = levels(resp), direction = ">", quiet = TRUE)
  if (as.numeric(pROC::auc(r2)) >= as.numeric(pROC::auc(r1))) r2 else r1
}

## ---- assemble scored, labelled rows ---------------------------------------
pr <- predict.lmvt(fitD, newdata = dat_fit, threshold = 0, innov_g = TRUE, innov_t = FALSE)
pr$cs   <- as.numeric(thr_s[as.character(pr$s)])
pr$risk <- 1 - pnorm((pr$cs - pr$muhat) / pr$sigmahat)          # exceedance risk score

df <- pr %>%
  inner_join(label_tbl, by = c("s","t")) %>%
  inner_join(dplyr::select(dat_fit, s, t, y), by = c("s","t")) %>%     # y = raw HR_mean
  filter(is.finite(risk)) %>%
  mutate(grp = factor(ifelse(label == 2, "stress", "non-stress"),
                      levels = c("non-stress","stress")))

## ---- headline: variance-aware risk score vs raw HR_mean -------------------
auc_risk <- as.numeric(pROC::auc(roc_best(df$grp, df$risk)))   # variance-aware score
auc_hr   <- as.numeric(pROC::auc(roc_best(df$grp, df$y)))      # raw HR_mean baseline

out <- data.frame(
  metric = c("Variance-aware risk score (Model D)", "Raw HR_mean baseline", "Difference"),
  AUC    = round(c(auc_risk, auc_hr, auc_risk - auc_hr), 3)
)
cat("Main WESAD result (stress vs non-stress; subject-specific c_s(0.70)):\n")
print(out, row.names = FALSE)
cat(sprintf("\nWindows used: %d (%d stress, %d non-stress) from %d participants.\n",
            nrow(df), sum(df$grp=="stress"), sum(df$grp=="non-stress"), length(unique(df$s))))

write.csv(out, file.path(OUTD, "wesad_main_risk_vs_hr.csv"), row.names = FALSE)
write.csv(df[, c("s","t","muhat","sigmahat","cs","risk","y","label")],
          file.path(OUTD, "wesad_main_scores.csv"), row.names = FALSE)
cat("\nWrote wesad_main_risk_vs_hr.csv and wesad_main_scores.csv to", OUTD, "\n")
cat("(Participant-clustered bootstrap CIs for these AUCs are computed in python/wesad_reanalyses.py.)\n")
