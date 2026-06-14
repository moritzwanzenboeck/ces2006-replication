# CES 2006 Replication

> _This README is auto-generated and maintained by Claude (Anthropic)._

Replication of **Cappiello, Engle & Sheppard (2006)**: *Asymmetric Dynamics in the Correlations of Global Equity and Bond Returns*, Journal of Financial Econometrics.

Extended sample period 2000–2025 (original: 1987–2002), 21 equity + 18 bond markets.

## Repository structure

| File / Folder | Description |
|---|---|
| `agdcc.R` | Three-stage QMLE estimator: DCC / ADCC / GDCC / AGDCC |
| `agdcc_core.cpp` | Rcpp/RcppArmadillo C++ kernels (Qt recursion, analytic scores, VCV) |
| `test_agdcc.R` | 120-assertion test suite |
| `0_data_msci_imi.R` | Parse MSCI IMI equity index files |
| `0_data_tr_gb.R` | Parse TR 10-year government bond index files |
| `1_data_manipulation.R` | FX conversion to USD |
| `2_returns_compute.R` | Thursday weekly log returns (T=1,355 × k=39) |
| `3_stage1_garch_fit.R` | BIC-selected univariate GARCH (8 specifications × 39 series) |
| `4_nonparametric_tests.R` | Variance and covariance asymmetry tests (NW-HAC) |
| `5_bootstrap.R` | Stationary bootstrap for median correlation significance |
| `6_dcc_estimation.R` | DCC/ADCC/GDCC/AGDCC estimation + 3-stage sandwich VCV |
| `7_tables.R` | All tables → `output/tables/*.tex` (incl. Table A1 data sources from `data_sources.xlsx`) |
| `8_charts.R` | All figures → `output/charts/*.pdf` |
| `data_sources.xlsx` | Index identifiers / provenance (Datastream + MSCI sheets); input to Table A1 |
| `mfe-toolbox/` | Sheppard's MFE MATLAB reference (`dcc.m`, `dcc_likelihood.m`, …) |
| `data/` | Processed data (weekly returns, GARCH residuals, test results) |
| `output/tables/` | LaTeX tables |
| `output/charts/` | LaTeX figures |
| `REFERENCE.md` | Comprehensive technical reference for all functions and design decisions |

## Data

- **Equity**: 21 MSCI IMI country indices (total return, gross), USD
- **Bond**: 18 Thomson Reuters 10-year Government Benchmark TR indices, USD
- **Returns**: Thursday-to-Thursday weekly log returns, 2000-01-03 to 2025-12-25
- `data/dcc_fits.RData` (234 MB) is not included in this repository due to size constraints

## Models

| Model | Parameters | Description |
|---|---|---|
| DCC | 2 | Scalar dynamic conditional correlation |
| ADCC | 3 | Scalar asymmetric DCC (asymmetric response to joint negative shocks) |
| GDCC | 2k=78 | Diagonal DCC (asset-specific parameters) |
| AGDCC | 3k=117 | Diagonal asymmetric DCC |

All models use full 3-stage sandwich standard errors (Engle & Sheppard 2001).

## Requirements

R packages: `rugarch`, `rmgarch`, `Rcpp`, `RcppArmadillo`, `numDeriv`, `sandwich`, `parallel`, `ggplot2`, `patchwork`, `tikzDevice`, `readxl`

Compile C++ kernels on first use: `Rcpp::sourceCpp("agdcc_core.cpp")` (inside `agdcc.R`).
