# 3_stage1_garch_specs.R
# Defines the 8 rugarch uGARCHspec objects used for Stage 1 BIC selection.
# AGARCH is excluded (not implemented in rugarch). CES 2006 Table 5 models:
#
#   CES 2006       rugarch model        rugarch submodel
#   GARCH          sGARCH               —
#   AVGARCH        fGARCH               AVGARCH
#   NARCH          fGARCH               NGARCH   (rugarch name for Higgins-Bera nonlinear ARCH)
#   EGARCH         eGARCH               —
#   ZARCH          fGARCH               TGARCH
#   GJR-GARCH      gjrGARCH             —
#   APARCH         apARCH               —
#   NAGARCH        fGARCH               NAGARCH
#
# Parameterization notes for Table 5 translation:
#   eGARCH:   alpha1 = sign effect (α), gamma1 = size effect (γ)
#   gjrGARCH: gamma1 = leverage multiplier; persistence = alpha1 + beta1 + 0.5*gamma1
#   apARCH:   gamma1 ∈ (-1,1) is leverage, delta is the power parameter
#   TGARCH:   eta11 is the asymmetry (Zakoian ZARCH)
#   NAGARCH:  eta21 is the shift (= -γ in Engle-Ng notation)
#
# Source this file to get a named list `garch_specs` (length 8).

packages <- c("rugarch")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

mean_spec   <- list(armaOrder = c(0, 0), include.mean = TRUE)
dist_spec   <- "norm"

garch_specs <- list(
  sGARCH = ugarchspec(
    variance.model = list(model = "sGARCH",   garchOrder = c(1, 1)),
    mean.model     = mean_spec,
    distribution.model = dist_spec
  ),
  AVGARCH = ugarchspec(
    variance.model = list(model = "fGARCH", submodel = "AVGARCH", garchOrder = c(1, 1)),
    mean.model     = mean_spec,
    distribution.model = dist_spec
  ),
  NGARCH = ugarchspec(
    variance.model = list(model = "fGARCH", submodel = "NGARCH",  garchOrder = c(1, 1)),
    mean.model     = mean_spec,
    distribution.model = dist_spec
  ),
  eGARCH = ugarchspec(
    variance.model = list(model = "eGARCH",  garchOrder = c(1, 1)),
    mean.model     = mean_spec,
    distribution.model = dist_spec
  ),
  TGARCH = ugarchspec(
    variance.model = list(model = "fGARCH", submodel = "TGARCH",  garchOrder = c(1, 1)),
    mean.model     = mean_spec,
    distribution.model = dist_spec
  ),
  gjrGARCH = ugarchspec(
    variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
    mean.model     = mean_spec,
    distribution.model = dist_spec
  ),
  apARCH = ugarchspec(
    variance.model = list(model = "apARCH",   garchOrder = c(1, 1)),
    mean.model     = mean_spec,
    distribution.model = dist_spec
  ),
  NAGARCH = ugarchspec(
    variance.model = list(model = "fGARCH", submodel = "NAGARCH", garchOrder = c(1, 1)),
    mean.model     = mean_spec,
    distribution.model = dist_spec
  )
)

cat("Loaded", length(garch_specs), "GARCH specs:", paste(names(garch_specs), collapse = ", "), "\n")
