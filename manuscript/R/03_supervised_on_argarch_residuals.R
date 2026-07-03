# ============================================================================
# 03_supervised_on_argarch_residuals.R
# WESAD: supervised classifiers (logistic, random forest, GBM) under
# leave-one-subject-out CV, on the window features alone vs. the features plus
# the AR-GARCH standardized residual.
# Run from the repository root:  Rscript manuscript/R/03_supervised_on_argarch_residuals.R
# ============================================================================

ROOT <- normalizePath(".")
EST  <- file.path(ROOT, "R", "temporalVarGuid.r")
DAT  <- file.path(ROOT, "manuscript", "data", "dat_fit.csv")
LAB  <- file.path(ROOT, "manuscript", "data", "label_tbl.csv")
OUTD <- file.path(ROOT, "manuscript", "outputs")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

source(EST)
suppressMessages({library(dplyr); library(pROC); library(randomForest); library(gbm)})

dat_fit   <- read.csv(DAT, check.names = FALSE)
label_tbl <- read.csv(LAB, check.names = FALSE)

## ---- 1) panel AR-GARCH on HR (no covariates) -> standardized residual -----
## Fit a covariate-free AR mean + GARCH variance via lmvt without a singular
## design: keep the real covariates but shrink ALL of them to zero with a large
## mean penalty (the AR lag is unpenalized), and exclude X from the variance.
suppressMessages(requireNamespace("glmnet"))
fit_ag <- lmvt(dat_fit, p = 1, q = 0, r = 1, s_ord = 1,
               use_x_in_variance = FALSE, lambda_beta = 1e3, standardize_X = TRUE)
pr_ag  <- predict.lmvt(fit_ag, newdata = dat_fit, threshold = 0, innov_g = TRUE)
z_df <- merge(pr_ag[, c("s","t","muhat","sigmahat")], dat_fit[, c("s","t","y")], by = c("s","t"))
z_df$z_argarch <- (z_df$y - z_df$muhat) / z_df$sigmahat
resid_tbl <- z_df[, c("s","t","z_argarch")]

## ---- 2) assemble supervised design ---------------------------------------
df <- dat_fit %>% inner_join(label_tbl, by = c("s","t")) %>%
  inner_join(resid_tbl, by = c("s","t")) %>%
  filter(label %in% c(0,1,2,3)) %>%
  mutate(grp = factor(ifelse(label == 2, "stress","non-stress"),
                      levels = c("non-stress","stress")))
feat_cols <- setdiff(names(dat_fit), c("s","t","y"))       # 77 features
X_base <- df[, feat_cols, drop = FALSE]
X_aug  <- cbind(X_base, z_argarch = df$z_argarch)

roc_best <- function(resp, score) {
  r1 <- pROC::roc(resp, score, levels = levels(resp), direction = "<", quiet = TRUE)
  r2 <- pROC::roc(resp, score, levels = levels(resp), direction = ">", quiet = TRUE)
  if (as.numeric(pROC::auc(r2)) >= as.numeric(pROC::auc(r1))) r2 else r1
}

loso_auc <- function(X) {
  subs <- sort(unique(df$s)); n <- nrow(df)
  ph <- data.frame(glm = rep(NA_real_, n), rf = NA_real_, gbm = NA_real_)
  set.seed(123)
  for (sid in subs) {
    tr <- df$s != sid; te <- df$s == sid
    yb <- as.integer(df$grp[tr] == "stress")
    ## logistic
    m1 <- suppressWarnings(glm(yb ~ ., data = data.frame(yb, X[tr,,drop=FALSE]), family = binomial()))
    ph$glm[te] <- as.numeric(predict(m1, newdata = X[te,,drop=FALSE], type = "response"))
    ## random forest
    m2 <- randomForest::randomForest(x = X[tr,,drop=FALSE], y = df$grp[tr],
                                     ntree = 500, mtry = floor(sqrt(ncol(X))))
    ph$rf[te] <- predict(m2, newdata = X[te,,drop=FALSE], type = "prob")[,"stress"]
    ## gbm
    m3 <- gbm::gbm(yb ~ ., data = data.frame(yb, X[tr,,drop=FALSE]),
                   distribution = "bernoulli", n.trees = 1500,
                   interaction.depth = 3, shrinkage = 0.01, n.minobsinnode = 10, verbose = FALSE)
    bi <- gbm::gbm.perf(m3, method = "OOB", plot.it = FALSE)
    ph$gbm[te] <- gbm::predict.gbm(m3, newdata = X[te,,drop=FALSE], n.trees = bi, type = "response")
  }
  sapply(ph, function(p) as.numeric(pROC::auc(roc_best(df$grp, p))))
}

auc_base <- loso_auc(X_base)
auc_aug  <- loso_auc(X_aug)
out <- data.frame(classifier = c("Logistic","Random forest","GBM"),
                  AUC_features_only      = round(auc_base,3),
                  AUC_features_plus_resid= round(auc_aug,3))
print(out)
write.csv(out, file.path(OUTD, "supervised_on_argarch_residuals.csv"), row.names = FALSE)
cat("\nWrote", file.path(OUTD, "supervised_on_argarch_residuals.csv"), "\n")
