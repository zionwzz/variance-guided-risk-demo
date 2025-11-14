# Variance-Aware Penalized Panel Modeling for Personalized Risk Assessment from Wearable Sensor Data

Minimal, reproducible code for a variance-guided penalized model with autoregressive dynamics for panel data. This repository provides tools for individualized risk assessment using time-series panel models with both mean and variance modeling.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Quick Start](#quick-start)
- [Model Specification](#model-specification)
- [Usage Example](#usage-example)
- [Outputs](#outputs)
- [Requirements](#requirements)
- [Code & Data Availability](#code--data-availability)
- [License](#license)
- [Citation](#citation)

---

## Overview

This repository implements a panel data model that simultaneously estimates:
- **Conditional mean** with autoregressive dynamics and distributed-lag effects
- **Conditional variance** with ARCH/GARCH dynamics and optional covariate effects
- **Individualized risk assessment** for exceedance probabilities

The model is particularly suited for wearable sensor data and other high-frequency panel applications.

---

## Repository Structure

```
.
├── R/
│   ├── simulator_scenarios.r      # Scenario generators (1–4)
│   └── temporalVarGuid1.r         # Core model: lmvt(), predict.lmvt()
├── scripts/
│   ├── 00_setup_env.R             # Seeds + single-thread settings + session info
│   ├── demo_algorithm.R           # One scenario, one fit (fast demo)
│   └── run_main_results.R         # Replication of paper-style figures
├── outputs/
│   └── .gitkeep                   # Figures/CSVs written here
├── LICENSE
├── .gitignore
└── README.md
```

---

## Quick Start

### 1. Install Required Packages

```r
install.packages(c("ggplot2", "dplyr", "tidyr"))
```

### 2. Run the Fast Demo (~1–2 minutes)

```bash
Rscript scripts/demo_algorithm.R
```

This demo:
- Simulates one scenario via `simulate_scenario_smallT()`
- Fits the variance-guided model with `use_x_in_variance = TRUE`
- Predicts $\hat{\mu}$ and $\hat{\sigma}$
- Computes exceedance risk at a pooled 90th-percentile cutoff
- Writes two PNG files to `outputs/`

### 3. Reproduce Main Simulation Results (longer runtime)

```bash
Rscript scripts/run_main_results.R
```

This script:
- Regenerates RMSE/Bias boxplots for $T \in \{120, 240\}$
- Evaluates across multiple scenarios and models A–D
- Saves raw and summary CSV files

---

## Model Specification

### Model Orders

The `lmvt()` function specifies the panel model using four order parameters **(p, q, r, s_ord)**:

| Parameter | Description |
|-----------|-------------|
| **p** | Autoregressive order in the **mean** (lags of $y$) |
| **q** | Distributed-lag order for exogenous predictors $x$ (in mean and variance) |
| **r** | ARCH order in the **variance** (lags of squared residuals $e^2$) |
| **s_ord** | GARCH order in the **variance** (lags of $\sigma^2$) |

**Note:** When `use_x_in_variance = TRUE`, the same lag order `q` is used for covariates in both the mean and variance equations.

### Model Decomposition

Let $y_{st}$ be the outcome for subject $s$ at time $t$. The model decomposes as:

$$y_{st} = \mu_{st} + \sigma_{st}\,\varepsilon_{st}$$

where:
- $\mathbb{E}(\varepsilon_{st} \mid x_{st}, \mathcal{F}_{t-1}) = 0$
- $\mathrm{Var}(\varepsilon_{st} \mid x_{st}, \mathcal{F}_{t-1}) = 1$ (standardized innovations)
- $\mathcal{F}_{t-1}$ is the information set up to time $t-1$
- $\mu_{st}$ is the conditional mean
- $\sigma_{st}$ is the conditional standard deviation

### Mean Equation (ARX: orders p, q)

The conditional mean with subject-specific intercepts $\alpha_s$, shared AR coefficients $\theta_i$, and distributed-lag effects $\beta_j$:

$$\mu_{st} = \alpha_s + \sum_{i=1}^{p} \theta_i\, y_{s,t-i} + \sum_{j=0}^{q} x_{s,t-j}^\top \beta_j$$

where:
- $x_{s,t-j}$ is the vector of exogenous predictors for subject $s$ at time $t-j$
- Coefficients $\{\beta_j\}$ are estimated under an $\ell_1$-penalty controlled by `lambda_beta`
- The response follows: $y_{st} = \mu_{st} + \sigma_{st}\varepsilon_{st}$

### Variance Equation (GARCH-X: orders s_ord, r, q)

Define the residuals $e_{st} := y_{st} - \mu_{st} = \sigma_{st}\varepsilon_{st}$. The **pure GARCH(s_ord, r)** component is:

```math
\delta_{st} := \omega_s + \sum_{\ell=1}^{s_{\text{ord}}} b_\ell\, \sigma_{s,t-\ell}^2 + \sum_{r'=1}^{r} a_{r'}\, e_{s,t-r'}^2
```

with constraints $\omega_s > 0$, $a_{r'} \geq 0$, $b_\ell \geq 0$ for positivity.

The **GARCH-X variance** augments this with lagged covariates:

```math
\sigma_{st}^2 = \delta_{st} + \sum_{m=0}^{q} x_{s,t-m}^\top \gamma_m
```

where:
- $\omega_s$ is a subject-specific baseline variance
- $\{b_\ell\}$ are GARCH parameters (effects of lagged variance $\sigma^2$)
- $\{a_{r'}\}$ are ARCH parameters (effects of squared residuals $e^2$)
- $\{\gamma_m\}$ are covariate effects in the variance, penalized by `lambda_gamma`

**Implementation note:** When `use_x_in_variance = TRUE`, covariates enter the variance equation with the same lag structure as the mean equation. To ensure $\sigma_{st}^2 > 0$, the implementation enforces non-negativity constraints and projects the linear predictor to a positive set when needed.

### Parameter Structure: Shared vs. Subject-Specific

**Subject-specific parameters:**
- $\alpha_s$ — intercept in the mean equation
- $\omega_s$ — baseline variance

**Shared population parameters:**
- $\{\theta_i\}$ — autoregressive coefficients
- $\{\beta_j\}$ — covariate effects in mean
- $\{\gamma_m\}$ — covariate effects in variance
- $\{a_{r'}\}$ — ARCH coefficients
- $\{b_\ell\}$ — GARCH coefficients

### Design Philosophy

**Shared lag structure for covariates:** When `use_x_in_variance = TRUE`, the same lag block of exogenous predictors $(x_{st}, x_{s,t-1}, \ldots, x_{s,t-q})$ appears in both mean and variance equations. This design choice:
- Reduces the number of tuning parameters
- Provides a parsimonious parameterization
- Improves estimation stability when $T$ is modest
- Aligns well with wearable-sensor panel settings where temporal dynamics are similar across equations

---

## Usage Example

### Prepare Your Data

Your data frame should include:
- `s` — subject ID (integer or factor)
- `t` — time index (integer)
- `y` — outcome variable (numeric)
- Additional columns — treated as exogenous predictors

### Fit the Model

```r
# Load the model code
source("R/temporalVarGuid1.r")

# Example: AR(1) mean + GARCH(1,1) variance with covariates
fit <- lmvt(
  df,
  p = 1,                    # AR order in mean
  q = 0,                    # Distributed-lag order for X (in both mean and variance)
  r = 1,                    # ARCH order in variance
  s_ord = 1,                # GARCH order in variance
  lambda_beta  = 0.01,      # Penalty for mean coefficients {beta_j}
  lambda_gamma = 0.01,      # Penalty for variance coefficients {gamma_m}
  standardize_X = TRUE,     # Standardize predictors
  use_x_in_variance = TRUE  # Include X in variance equation
)
```

### Make Predictions

```r
# Obtain fitted values and innovations
pred <- predict.lmvt(fit, newdata = df, innov_g = TRUE)

# Components:
# pred$muhat      - fitted conditional mean
# pred$sigmahat   - fitted conditional standard deviation
# pred$innov      - standardized innovations
```

### Calculate Exceedance Risk

```r
# Define a risk threshold (e.g., 90th percentile of observed y)
c_cut <- quantile(df$y, 0.90, na.rm = TRUE)

# Compute exceedance probability P(y > c_cut | F_{t-1})
z <- (c_cut - pred$muhat) / pmax(pred$sigmahat, 1e-8)
p_exceed <- 1 - pnorm(z)
```

---

## Outputs

### Demo Outputs (subject-level visualizations)
- `outputs/scen<id>_subj1_y_mu_band.png` — Time series with fitted mean and uncertainty bands
- `outputs/scen<id>_subj1_risk_vs_hits.png` — Predicted risk vs. actual exceedances

### Main Results (simulation study)
Boxplots comparing models A–D:
- `outputs/ABCD_s15_T120_ar1_RMSE.png`
- `outputs/ABCD_s15_T120_ar1_Bias.png`
- `outputs/ABCD_s15_T240_ar1_RMSE.png`
- `outputs/ABCD_s15_T240_ar1_Bias.png`

### CSV Files
Raw and summarized simulation results:
- `outputs/rep_raw_s15t120_by_rep_30reps.csv`
- `outputs/rep_raw_s15t240_by_rep_30reps.csv`
- `outputs/rep_summary_s15t120_30reps.csv`
- `outputs/rep_summary_s15t240_30reps.csv`

---

## Requirements

### Software
- **R ≥ 4.2**

### R Packages
- **ggplot2** — visualization
- **dplyr** — data manipulation
- **tidyr** — data tidying

### Reproducibility Setup

The script `scripts/00_setup_env.R`:
- Fixes a global RNG seed for reproducibility
- Caps BLAS/LAPACK to single thread for determinism
- Writes `session_info.txt` in the repository root

### Optional: Lock Dependencies with renv

For enhanced reproducibility:

```r
install.packages("renv")
renv::init()
renv::snapshot()  # Creates renv.lock for reproducible package restoration
```

---

## Code & Data Availability

**Code:** All code to implement the method and reproduce manuscript figures is included in this repository.
- Minimal example: `scripts/demo_algorithm.R`
- Full replication: `scripts/run_main_results.R`

**Data:** All analyses use simulated data generated by `R/simulator_scenarios.r`. No external datasets are required.

---

## License

Released under the MIT License. See `LICENSE` for details.

---

## Citation

If you use this repository in your research, please cite the accompanying manuscript.

**Note:** A DOI (Zenodo/OSF) and `CITATION.cff` will be added upon release.

---

## Contact & Support

For questions or issues, please open an issue on the GitHub repository.
