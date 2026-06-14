# Project Reference: CES 2006 Replication

**Last updated:** 2026-06-13  
**Author:** Moritz Wanzenböck (UZH, FS 2026)

> _This reference document is auto-generated and maintained by Claude (Anthropic)._

Replication of Cappiello, Engle & Sheppard (2006, *JFEC*) over 2000–2025.  
Implements the full empirical pipeline: GARCH selection → asymmetry tests → bootstrap → DCC/ADCC/GDCC/AGDCC estimation → tables and charts.

**Key adaptations vs CES 2006:**
- Period extended to 2025 (T = 1355 weekly observations)
- 18 bond markets (vs 13 in CES)
- 10-year TR indices only; no structural break (4 models, not 12)

---

## Repository Layout

```
ces2006-replication/
├── agdcc_core.cpp          # C++ kernels (Rcpp/RcppArmadillo)
├── agdcc.R                 # Three-stage estimator (main R implementation)
├── test_agdcc.R            # 120-assertion test suite
│
├── 0_data_msci_imi.R       # Ingest MSCI IMI Excel files (raw data not included)
├── 0_data_tr_gb.R          # Ingest TR 10-yr Government Benchmark CSVs (raw data not included)
├── 1_data_manipulation.R   # FX conversion to USD
├── 2_returns_compute.R     # Thu-to-Thu weekly log returns
├── 2_returns_inspect.R     # Diagnostic plots (run manually)
├── 3_stage1_garch_specs.R  # Define 8 rugarch spec objects
├── 3_stage1_garch_fit.R    # Fit 8 specs × 39 series; BIC selection
├── 4_nonparametric_tests.R # Variance + covariance asymmetry tests
├── 5_bootstrap.R           # Stationary bootstrap (median correlation)
├── 6_dcc_estimation.R      # Fit four DCC models; compute 3-stage VCV
├── 7_tables.R              # All LaTeX tables → output/tables/
├── 8_charts.R              # All PDF figures  → output/charts/
│
├── data/                   # Processed data (weekly returns, GARCH fits, test results)
├── output/
│   ├── tables/             # .tex table files
│   └── charts/             # .pdf + TikZ .tex figure files
│
├── mfe-toolbox/            # Sheppard MATLAB reference (dcc.m, etc.)
├── README.md               # Repository overview
└── REFERENCE.md            # ← This file; comprehensive project reference
```

---

## Execution Order

Scripts 0a–2 require proprietary raw data files (MSCI IMI Excel, Thomson Reuters CSVs) not included in this repository. All downstream steps use the pre-computed `data/` files provided here.

| # | Script | Runtime | Output |
|---|--------|---------|--------|
| 0a | `0_data_msci_imi.R` | fast | `data/msci_*.RData` |
| 0b | `0_data_tr_gb.R` | fast | `data/tr_gb_*.RData` |
| 1 | `1_data_manipulation.R` | fast | `data/prices_*.RData` |
| 2 | `2_returns_compute.R` | fast | `data/weekly_returns_usd.RData` |
| 3 | `3_stage1_garch_fit.R` | ~10 min | `data/stage1_fits.RData` |
| 4 | `4_nonparametric_tests.R` | ~2 min | `data/nonparam_tests.RData` |
| 5 | `5_bootstrap.R` | ~5 min | `data/bootstrap_corr.RData` |
| **6** | **`6_dcc_estimation.R`** | **~11 h** | **`data/dcc_fit_*_3stage_tmp.rds`** |
| 7 | `7_tables.R` | fast | `output/tables/*.tex` |
| 8 | `8_charts.R` | ~1 min | `output/charts/*.{pdf,tex}` |

Steps 3–8 can be run directly using the data files provided in this repository. Step 6 can be skipped — the `dcc_fit_*_3stage_tmp.rds` output files are included.

---

## Data

### Assets (k = 39)

- **Equities (21):** MSCI IMI country TR indices; prefix `eq_`
  - Australasia: `eq_australia`, `eq_hong_kong`, `eq_japan`, `eq_new_zealand`, `eq_singapore`
  - Europe (13): `eq_austria`, `eq_belgium`, `eq_denmark`, `eq_france`, `eq_germany`, `eq_ireland`, `eq_italy`, `eq_netherlands`, `eq_norway`, `eq_spain`, `eq_sweden`, `eq_switzerland`, `eq_united_kingdom`
  - N. America: `eq_canada`, `eq_mexico`, `eq_usa`
- **Bonds (18):** Thomson Reuters 10-yr Government Benchmark TR indices; prefix `bd_`
  - Countries: Australia, Austria, Belgium, Canada, Denmark, France, Germany, Ireland, Italy, Japan, Netherlands, New Zealand, Norway, Spain, Sweden, Switzerland, United Kingdom, United States

### Data files

| File | Contents |
|------|----------|
| `data/weekly_returns_usd.RData` | T×k matrix `returns_mat` of Thu-to-Thu USD log returns |
| `data/stage1_fits.RData` | `stage1_best`, `stage1_selected`, `e_std` (T×k), `h_mat` (T×k) |
| `data/nonparam_tests.RData` | `var_results` (39×4), `test_results` (741×6), `table4_summary` |
| `data/bootstrap_corr.RData` | `obs_medians`, `boot_medians` (B×3), `p_values`, `params` |
| `data/dcc_fit_{MODEL}_none_tmp.rds` | Phase-1 fit (optimised parameters, no VCV) for DCC/ADCC/GDCC/AGDCC |
| `data/dcc_fit_{MODEL}_3stage_tmp.rds` | Phase-2 fit with full 3-stage VCV and SEs for DCC/ADCC/GDCC/AGDCC |

`data/dcc_fits.RData` (all four models combined, 234 MB) is not included due to size. It can be reconstructed from the `*_3stage_tmp.rds` files via the bundling step at the end of `6_dcc_estimation.R`.

---

## `agdcc_core.cpp` — C++ Kernels

Compiled once per session via `Rcpp::sourceCpp("agdcc_core.cpp")`.  
Uses **RcppArmadillo** for linear algebra and **OpenMP** for within-call parallelism.

### Model integer codes

```cpp
1 = DCC_SCALAR   2 = ADCC_SCALAR   3 = GDCC_DIAG   4 = AGDCC_DIAG
```

### Exported functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `qt_recursion_cpp` | `(psi_D, P, N, ee, nn, bc, bca, model) → List(Qt, Rt)` | Full Qt/Rt recursion; returns k×k×T cubes |
| `stage3_qll_fast_cpp` | `(psi_D, P, N, e_std, ee, nn, bc, bca, model) → double` | Stage-3 quasi-LL in one pass, no cube allocation |
| `analytic_scores_parallel_cpp` | `(psi_D, P, N, e_std, ee, nn, bc, bca, model, Qt, Rt) → List(scores, B22)` | Analytic Ding-Engle score recursion, OpenMP-parallelised over np parameters |
| `ewma_backcast_cpp` | `(X, lambda=0.94) → mat` | EWMA backcast of a k×k×T cube → k×k starting value |
| `stage3_qll_grid_cpp` | `(grid, P, N, e_std, ee, nn, bc, bca, model, n_cores) → vec` | Evaluate LL for each row of starting-value grid |
| `build_outer_cubes_cpp` | `(e_std, n_std) → List(ee, nn)` | Build k×k×T outer-product cubes from standardised residuals |
| `compute_A32_cpp` | `(psi_D, P, N, e_std, ee, nn, bc, bca, model, n_cores) → mat` | A₃₂ cross-Jacobian via OpenMP over n₂ pairs (n₃×n₂ result) |
| `qt_and_mean_scores_cpp` | `(psi_D, P, N, e_std, ee, nn, bc, bca, model, n_cores) → vec` | Qt recursion + mean scores (1/T·Σ ∂ℓ/∂ψ), np-vector; used for optimizer gradient and A₃₃ Jacobian |

### Qt recursion (all models)

```
Qt = (P̄ − A'P̄A − B'P̄B − G'N̄G) + A'(ε̂_{t-1}ε̂'_{t-1})A + G'(η_{t-1}η'_{t-1})G + B'Q_{t-1}B
Rt = diag(Qt)^{-1/2} Qt diag(Qt)^{-1/2}
ℓ_t = −½(log|Rt| + ε̂_t' Rt⁻¹ ε̂_t − ε̂_t'ε̂_t)
```

The intercept `P̄ − A'P̄A − B'P̄B − G'N̄G` ensures E[Qt] = P̄ at stationarity (variance targeting).

Scalar models: A=aI, G=gI, B=bI.  
Diagonal models: A=diag(a₁…aₖ), G=diag(g₁…gₖ), B=diag(b₁…bₖ).  
Parameters stored as raw values; A, G, B formed as `psi²` (squared) in the recursion.

---

## `agdcc.R` — Three-Stage Estimator

**Source once:** `source("agdcc.R")` (auto-compiles C++ if needed).

### Constants

```r
DCC_SCALAR = 1L;  ADCC_SCALAR = 2L
GDCC_DIAG  = 3L;  AGDCC_DIAG  = 4L
```

### Utility functions

| Function | Purpose |
|----------|---------|
| `covnw(X, L=NULL)` | Newey-West HAC covariance of T×p score matrix; bandwidth `min(floor(1.2*T^(1/3)), T)` (matches Sheppard's `covnw.m`) |
| `.pairs_lower(k)` | K(K-1)/2 × 2 matrix of (i>j) pairs, column-major (R̄ block ordering) |
| `.pairs_lower_diag(k)` | K(K+1)/2 × 2 matrix of (i≥j) pairs, column-major (N̄ block ordering) |
| `.safe_solve(A)` | Solve with ridge fallback (1e-8) then `MASS::ginv` |
| `.sym(M)` | Force numerical symmetry: `(M + t(M)) / 2` |
| `is_pd(M)` | TRUE if Cholesky succeeds |
| `compute_delta(P, N)` | Max eigenvalue of `P^{-1/2} N P^{-1/2}` (stationarity scalar for ADCC/AGDCC) |

### Stage 1: `fit_stage1(data, stage1_fit=NULL, n_cores=1L, return_fits=FALSE)`

Extracts per-obs conditional variances and standardised residuals from univariate GARCH fits.

- **`stage1_fit = NULL`:** fits internal sGARCH(1,1) per asset (for testing/standalone use)
- **`stage1_fit = ugarchMultifit`:** extracts from rugarch multifit object
- **`stage1_fit = list<ugarchfit>`:** extracts from pre-fitted list (production path from `3_stage1_garch_fit.R`)
- **`return_fits = TRUE`:** also stores `$specs` and `$thetas` — required for 3-stage VCV

**Returns:** `list(H, e_std, [specs, thetas])`

### Stage 2: `compute_intercepts(e_std)`

Computes variance-targeting intercepts:
- `P = cor(e_std)` — unconditional correlation matrix P̄
- `N = E[η_t η_t']` where `η_t = ε̂_t * (ε̂_t < 0)` — asymmetric intercept N̄

**Returns:** `list(P, N, n_std)`

### Stage 3: optimiser

| Function | Purpose |
|----------|---------|
| `stage3_grid_search(...)` | Multi-start grid: 9 candidates for DCC/GDCC, 18 for ADCC/AGDCC |
| `optimise_stage3(...)` | L-BFGS-B → Nelder-Mead (scalar, np≤5) → L-BFGS-B; analytical gradient |
| `stage3_obj(psi, ...)` | Objective = Stage-3 QL + stationarity penalty (activated at 0.9999) |

**Gradient:** `gr(ψ) = −T · qt_and_mean_scores_cpp(ψ, …) + penalty_grad`  
The T factor is critical — `qt_and_mean_scores_cpp` returns MEAN scores (1/T·Σ), but the objective is the SUM.

### VCV: Full 3-block sandwich (Sheppard 2012 MFE toolbox)

**Block structure:**
```
A₀ = [A₁₁   0    0  ]     n₁ = Σᵢ nᵢ (GARCH params, typically 5 per asset)
     [A₂₁  A₂₂   0  ]     n₂ = K(K-1)/2 [+ K(K+1)/2 for ADCC/AGDCC]
     [A₃₁  A₃₂  A₃₃ ]     n₃ = length(psi_D)

VCV(ψ̂) = [A₀⁻¹ B₀ (A₀⁻¹)' / T]_{(n₁+n₂+1):n_tot, (n₁+n₂+1):n_tot}
```

where B₀ = covnw([S₁ | S₂ | S₃]) is the Newey-West HAC covariance of the joint T×n_tot score matrix.

| Block | What | How |
|-------|------|-----|
| A₁₁ | GARCH info (block-diagonal) | Empirical information: `crossprod(S1[:,i]) / T` per asset |
| A₂₂ | Intercept pseudo-Hessian | Exactly `I_{n₂}` (analytic identity) |
| A₃₃ | DCC dynamics info | `numDeriv::jacobian(qt_and_mean_scores_cpp, psi_D)` |
| A₂₁ | Intercept ← GARCH cross | Forward finite differences on `ugarchfilter`; zero at solution (m_base=0) |
| A₃₁ | DCC ← GARCH cross | Forward FD on `ugarchfilter` + `qt_and_mean_scores_cpp` |
| A₃₂ | DCC ← intercept cross | `compute_A32_cpp` (C++ OpenMP over n₂ pairs) |
| S₁ | GARCH per-obs scores | Forward FD on per-obs GARCH log-likelihoods |
| S₂ | Intercept per-obs scores | Analytic: `ε̂ᵢε̂ⱼ/sᵢⱼ − P̄ᵢⱼ` (R̄ block); `ηᵢηⱼ − N̄ᵢⱼ` (N̄ block) |
| S₃ | DCC per-obs scores | `analytic_scores_parallel_cpp` (Ding-Engle recursion) |

**`vcv_3stage` return list:**
```r
list(VCV, se, A11, A21, A31, A32, A33, B0, n1, n2, n3, scores, method="3stage")
```

**Rank-deficiency:** when n_tot ≥ T, B₀ is rank-deficient. A warning is emitted but the ψ-block SEs remain valid (S₃ has rank n₃ ≪ T).

| Model | K=39 n_tot | Feasible? |
|-------|-----------|-----------|
| DCC   | 938       | ✓ (< T=1355) |
| ADCC  | 1719      | ✗ (rank-deficient, warns) |
| GDCC  | 1014      | ✓ |
| AGDCC | 1833      | ✗ (warns) |

### Main entry point: `fit_agdcc(data, model, stage1_fit, vcv_method, n_cores, cpp_file)`

```r
fit_agdcc(data,
  model      = "AGDCC",          # "DCC" | "ADCC" | "GDCC" | "AGDCC"
  stage1_fit = NULL,             # NULL | ugarchMultifit | list<ugarchfit>
  vcv_method = c("none","3stage"),
  n_cores    = 1L,
  cpp_file   = "agdcc_core.cpp")
```

**`agdcc_fit` object fields:**

| Field | Type | Contents |
|-------|------|----------|
| `psi_D` | named vec | DCC parameter estimates |
| `ll` | scalar | Stage-3 quasi-log-likelihood |
| `Ht` | k×k×T array | Full conditional covariance matrices |
| `Rt` | k×k×T array | Conditional correlation matrices |
| `Qt` | k×k×T array | Qt recursion output |
| `VCV` | n3×n3 mat | VCV of ψ̂ (NULL if vcv="none") |
| `se` | named vec | Standard errors (NA if vcv="none") |
| `tstat` | named vec | t-statistics |
| `pval` | named vec | Two-sided p-values |
| `H_univ` | T×k mat | Univariate GARCH conditional variances |
| `e_std` | T×k mat | GARCH standardised residuals |
| `P` | k×k mat | Unconditional correlation matrix P̄ |
| `N` | k×k mat | Asymmetric intercept matrix N̄ |
| `delta` | scalar | Stationarity scalar (max eigenvalue) |
| `model` | string | "DCC" / "ADCC" / "GDCC" / "AGDCC" |
| `vcv_method` | string | "none" / "3stage" |
| `T`, `k` | ints | Dimensions |
| `vcv_detail` | list | Full VCV blocks (A11, A21, A31, A32, A33, B0) |

### S3 methods

```r
print.agdcc_fit(x)        # parameter table with stars
AIC.agdcc_fit(x, k=2)     # AIC = -2*ll + k*np; also BIC (k=log(T))
coef.agdcc_fit(x)         # named psi_D vector
vcov.agdcc_fit(x)         # VCV matrix
```

### Parameter naming convention

- Scalar (DCC): `"a"`, `"b"`
- Scalar (ADCC): `"a"`, `"g"`, `"b"`
- Diagonal (GDCC): `"[asset_name].a1"`, `"[asset_name].b1"`
- Diagonal (AGDCC): `"[asset_name].a1"`, `"[asset_name].g1"`, `"[asset_name].b1"`

---

## `test_agdcc.R` — Test Suite

Run with `Rscript test_agdcc.R`. **120 assertions** across 13 test groups:

| Group | What |
|-------|------|
| T1 | C++ exports callable |
| T2 | Qt/Rt: symmetry, PD, unit diagonal; all 4 models |
| T3 | Stage-3 LL finite for valid params; 1e9 for violations |
| T4 | Analytic vs numerical scores: relative error < 10% |
| T5 | `fit_agdcc()` all 4 models with `vcv_method="none"` |
| T6 | `stage1_fit` modes: ugarchMultifit, list\<ugarchfit\> |
| T7 | Stationarity constraint satisfied at convergence |
| T8 | VCV: square, symmetric (tol 1e-8), PSD |
| T10 | S3 methods: print, AIC/BIC, coef, vcov |
| T11 | Conditional correlations ∈ (−1, 1) for all t, i, j |
| T12 | Serial == parallel scores at n_cores = 1 |
| T13 | 3-stage VCV: positive SEs, PSD VCV, A11/A21/A31/A32 present |

Test data: `set.seed(42); x <- matrix(rnorm(500*3), 500, 3)` (k=3, T=500).

---

## `6_dcc_estimation.R` — DCC Estimation Script

Fits all four models using a two-phase workflow:

**Phase 1** — optimisation only (output provided as `dcc_fit_*_none_tmp.rds`):
```r
fits[[m]] <- fit_agdcc(returns_mat, model=m, stage1_fit=stage1_best,
                        vcv_method="none", n_cores=...)
saveRDS(fits[[m]], sprintf("data/dcc_fit_%s_none_tmp.rds", m))
```

**Phase 2** — 3-stage VCV added post-hoc (output provided as `dcc_fit_*_3stage_tmp.rds`):
```r
fit <- readRDS(sprintf("data/dcc_fit_%s_none_tmp.rds", m))
fit <- compute_vcv(fit, vcv_method="3stage", n_cores=...)
saveRDS(fit, sprintf("data/dcc_fit_%s_3stage_tmp.rds", m))
```

### Full-sample results (T=1355, k=39)

| Model | np | LL | BIC |
|-------|----|----|-----|
| DCC   | 2  | 14779 | −29543 |
| ADCC  | 3  | 14862 | −29707 |
| GDCC  | 78 | 16143 | −31818 |
| AGDCC | 117 | 16279 | −31955 |

---

## `7_tables.R` — LaTeX Table Generation

Outputs to `output/tables/`. Tables 6b, 7, 8 require `data/dcc_fits.RData` (not provided). Load the included `data/dcc_fit_*_3stage_tmp.rds` files and bundle them into `dcc_fits.RData` using the last block of `6_dcc_estimation.R`.

### Tables

| File | Description | Requires |
|------|-------------|---------|
| `table1_descriptive.tex` | Annualised mean/SD, skew, kurtosis, JB test for all 39 series | stage1_fits |
| `table2_avg_correlations.tex` | Mean/min/max correlations + 3×3 regional matrices | weekly_returns |
| `table3a_corr_equity.tex` | 21×21 equity correlation matrix (lower-tri) | weekly_returns |
| `table3b_corr_eq_bond.tex` | 21×18 equity-bond cross-correlation matrix | weekly_returns |
| `table3c_corr_bond.tex` | 18×18 bond correlation matrix (lower-tri) | weekly_returns |
| `table4_partial_cov.tex` | Proportion of pairs with significant covariance asymmetry | nonparam_tests |
| `table5_garch_params.tex` | Stage-1 GARCH model selection and parameter estimates | stage1_fits |
| `table6a_dcc_params.tex` | Diagonal DCC parameter estimates (GDCC/AGDCC per-asset, DCC/ADCC scalar) | dcc_fits |
| `table6b_dcc_comparison.tex` | Model comparison: np, LL, BIC | dcc_fits |
| `table7_equity_var_corr.tex` | Equity vol–correlation relationship (AGDCC) | dcc_fits |
| `table8_bond_var_corr.tex` | Bond vol–correlation relationship (AGDCC) | dcc_fits |
| `table_a1_data_sources.tex` | Data sources / index identifiers: country, index name, Datastream mnemonic (bonds), provider code (investing.com / MSCI), in-CES-2006 flag (bonds), source URL. Two-panel (Equities/Bonds) `longtable` | `data_sources.xlsx` |

**Table 1 star convention (reversed — stars mean approximately normal):**
- `\text{*}` — p > 0.05 (fail to reject normality at 5%)
- `\text{**}` — 0.01 < p ≤ 0.05 (fail to reject at 1%)
- no star — p ≤ 0.01 (strongly rejects)

**Table 6a star convention:** `\text{*}` = insignificant at 5% (`|z| < 1.96`).

---

## `8_charts.R` — Figure Generation

Outputs to `output/charts/`. Figures 2b, 3, 6–10 require `data/dcc_fits.RData` (not provided; bundle from `dcc_fit_*_3stage_tmp.rds` first).

| File | Description |
|------|-------------|
| `fig1_volatility_nic.pdf/.tex` | Volatility NIC curves (2 panels: equities / bonds) |
| `fig2a_corr_nic_3d.pdf/.tex` | 3D correlation NIC surface (pgfplots) |
| `fig2b_corr_nic_views.pdf/.tex` | 2D NIC cross-sections (four views) |
| `fig3_covariance_nic_3d.pdf/.tex` | 3D covariance NIC surface (return units ±0.05) |
| `fig4_equity_volatility.pdf/.tex` | Equity conditional volatilities (AGDCC) |
| `fig5_bond_volatility.pdf/.tex` | Bond conditional volatilities; uses stage1_fits only |
| `fig6_equity_correlations.pdf/.tex` | Pairwise correlations: France/Germany/Italy/UK |
| `fig7_regional_equity_corr.pdf/.tex` | Regional equity index correlations |
| `fig8_bond_correlations.pdf/.tex` | Bond regional correlations |
| `fig9_bond_corr_selected.pdf/.tex` | Selected bond pair correlations |
| `fig10_equity_bond_corr.pdf/.tex` | Equity-bond correlations |

---

## Regional Groupings

```r
reg$eq_australasia  # Australia, Hong Kong, Japan, New Zealand, Singapore (5)
reg$eq_europe       # Austria, Belgium, Denmark, France, Germany, Ireland, Italy,
                    # Netherlands, Norway, Spain, Sweden, Switzerland, UK (13)
reg$eq_emu          # Austria, Belgium, France, Germany, Ireland, Italy,
                    # Netherlands, Spain (8)
reg$eq_namerica     # Canada, Mexico, USA (3)
reg$bd_australasia  # Australia, Japan, New Zealand (3)
reg$bd_europe       # Austria, Belgium, Denmark, France, Germany, Ireland, Italy,
                    # Netherlands, Norway, Spain, Sweden, Switzerland, UK (13)
reg$bd_namerica     # Canada, United States (2)
```

---

## MFE Toolbox Reference (`mfe-toolbox/`)

Sheppard's MATLAB reference implementation. Key files:

| File | Role |
|------|------|
| `dcc.m` | Main estimation; 3-stage inference path (lines 420–471) |
| `dcc_likelihood.m` | Stage-3 log-likelihood |
| `dcc_inference_objective.m` | Score functions for inference |
| `dcc_fit_variance.m` | Stage-1 variance fitting |
| `covnw.m` | NW-HAC; bandwidth `min(floor(1.2*T^(1/3)),T)` — matched in R `covnw()` |
| `gradient_2sided.m` | Central-difference gradient |
| `hessian_2sided_nrows.m` | Numerical Hessian (nrows-only variant) |

**Key difference from Sheppard:** this implementation uses analytic Stage-3 scores (`analytic_scores_parallel_cpp`) and analytic A₂₂=I; Sheppard uses numerical differentiation for both.
