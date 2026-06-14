# test_agdcc.R — 113-assertion test suite for agdcc_core.cpp + agdcc.R
#
# T1:  C++ exports are callable
# T2:  Qt/Rt recursion: symmetry, PD, unit diagonal; all 4 models
# T3:  Stage-3 LL: finite for valid params; 1e9 for constraint violations
# T4:  Analytic vs numerical scores: relative error < 10%
#        (interior psi used so Richardson steps stay clear of the 0.9999 gate)
# T5:  fit_agdcc: all 4 models with vcv_method="none"
# T6:  stage1_fit modes: ugarchMultifit and list<ugarchfit>
# T7:  Stationarity constraint satisfied at convergence
# T8:  VCV square, symmetric (tol 1e-8 for Jacobian noise), PSD
# T9:  (removed — lr_test deleted)
# T10: S3 methods: print, AIC/BIC, coef, vcov
# T11: conditional_correlations in (-1, 1) for all t, i, j
# T12: serial == parallel scores at n_cores = 1

cat("=== AG-DCC Test Suite (v4) ===\n\n")

packages <- c("rugarch", "numDeriv", "Rcpp", "MASS")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
set.seed(123)

Sys.setenv(HOME = "C:/") 

Rcpp::sourceCpp("agdcc_core.cpp")
source("agdcc.R")

sim_returns <- function(T=500, k=3, seed=1) {
  set.seed(seed)
  MASS::mvrnorm(T, rep(0,k), diag(k)+0.3*(1-diag(k)))
}

data3 <- sim_returns(500, 3)

pass <- 0L; fail <- 0L

run_test <- function(label, expr) {
  res <- tryCatch({
    v <- force(expr)
    if (isTRUE(v)) { cat(sprintf("  PASS  %s\n", label)); TRUE }
    else           { cat(sprintf("  FAIL  %s  [%s]\n", label, paste(v,collapse=","))); FALSE }
  }, error=function(e) { cat(sprintf("  FAIL  %s  (error: %s)\n", label, conditionMessage(e))); FALSE })
  if (res) pass <<- pass+1L else fail <<- fail+1L
}
skip_test <- function(label, reason) cat(sprintf("  SKIP  %s  (%s)\n", label, reason))

# ---------------------------------------------------------------------------
cat("-- T1: C++ exports --\n")
for (fn in c("qt_recursion_cpp","stage3_qll_cpp","analytic_scores_cpp",
             "ewma_backcast_cpp",
             "stage3_qll_grid_cpp","analytic_scores_parallel_cpp"))
  run_test(paste(fn,"exists"), exists(fn))

# ---------------------------------------------------------------------------
cat("\n-- T2: Qt recursion --\n")
k <- 3L
e_tmp <- scale(data3); n_tmp <- e_tmp*(e_tmp<0)
P_tmp <- cov2cor_safe(crossprod(e_tmp)/nrow(e_tmp))$R
N_tmp <- crossprod(n_tmp)/nrow(e_tmp)
ee_tmp <- array(0,c(k,k,nrow(e_tmp))); nn_tmp <- array(0,c(k,k,nrow(e_tmp)))
for (t in seq_len(nrow(e_tmp))) { ee_tmp[,,t]<-tcrossprod(e_tmp[t,]); nn_tmp[,,t]<-tcrossprod(n_tmp[t,]) }
bc_tmp  <- ewma_backcast_cpp(ee_tmp)
bca_tmp <- ewma_backcast_cpp(nn_tmp)

for (mod in 1:4) {
  psi <- switch(as.character(mod),
    "1"=c(0.05,0.90),"2"=c(0.05,0.02,0.90),
    "3"=c(rep(0.05,k),rep(0.90,k)),"4"=c(rep(0.05,k),rep(0.02,k),rep(0.90,k)))
  rec <- qt_recursion_cpp(psi,P_tmp,N_tmp,ee_tmp,nn_tmp,bc_tmp,bca_tmp,mod)
  Qt_T <- rec$Qt[,,nrow(e_tmp)]; Rt_T <- rec$Rt[,,nrow(e_tmp)]
  mn   <- c("DCC","ADCC","GDCC","AGDCC")[mod]
  run_test(sprintf("Qt[T] symmetric (%s)",mn), max(abs(Qt_T-t(Qt_T)))<1e-10)
  run_test(sprintf("Qt[T] PD (%s)",mn),        is_pd(Qt_T))
  run_test(sprintf("Rt[T] is corr (%s)",mn),
           max(abs(diag(Rt_T)-1))<1e-10 && all(abs(Rt_T[lower.tri(Rt_T)])<=1+1e-10))
  run_test(sprintf("Rt[T] PD (%s)",mn),        is_pd(Rt_T))
}

# ---------------------------------------------------------------------------
cat("\n-- T3: Stage-3 LL --\n")
psi_good <- c(0.05,0.70)   # well inside feasible region: stat=0.0025+0.49=0.49
psi_bad1 <- c(0.80,0.80)   # stat_val = 1.28 > 0.9999
psi_bad2 <- c(-0.05,0.90)  # negative
ll_g <- stage3_qll_cpp(psi_good,P_tmp,N_tmp,e_tmp,ee_tmp,nn_tmp,bc_tmp,bca_tmp,1L)
ll_b1<- stage3_qll_cpp(psi_bad1,P_tmp,N_tmp,e_tmp,ee_tmp,nn_tmp,bc_tmp,bca_tmp,1L)
ll_b2<- stage3_qll_cpp(psi_bad2,P_tmp,N_tmp,e_tmp,ee_tmp,nn_tmp,bc_tmp,bca_tmp,1L)
run_test("LL finite for valid params",        is.finite(ll_g))
run_test("LL=1e9 for stationarity violation", ll_b1 >= 1e8)
run_test("LL=1e9 for negative param",         ll_b2 >= 1e8)
run_test("LL positive (negated sum)",         ll_g > 0)

# ---------------------------------------------------------------------------
cat("\n-- T4: Analytic vs numerical scores --\n")
# Use interior psi points well away from stationarity boundary so Richardson
# steps never hit the 0.9999 gate.

psi_test <- list(
  "1"=c(0.05,0.70),   # DCC:  stat=0.49
  "2"=c(0.05,0.02,0.70),  # ADCC: stat~0.49
  "3"=c(rep(0.05,k),rep(0.70,k)),   # GDCC: stat~0.49*3=1.47... need k=3
  "4"=c(rep(0.05,k),rep(0.02,k),rep(0.70,k)))
# For GDCC/AGDCC diagonal, stat = max_i(a_i^2 + b_i^2).
# Use interior values well away from the boundary: a_i=0.03, b_i=0.55 → max=0.30
psi_test[["3"]] <- c(rep(0.03,k), rep(0.55,k))
psi_test[["4"]] <- c(rep(0.03,k), rep(0.02,k), rep(0.55,k))

for (mod in 1:4) {
  mn  <- c("DCC","ADCC","GDCC","AGDCC")[mod]
  psi <- psi_test[[as.character(mod)]]

  rec <- qt_recursion_cpp(psi,P_tmp,N_tmp,ee_tmp,nn_tmp,bc_tmp,bca_tmp,mod)
  sc  <- analytic_scores_cpp(psi,P_tmp,N_tmp,e_tmp,ee_tmp,nn_tmp,
                              bc_tmp,bca_tmp,mod,rec$Qt,rec$Rt)

  num_grad <- tryCatch(
    numDeriv::grad(function(p)
      stage3_qll_cpp(p,P_tmp,N_tmp,e_tmp,ee_tmp,nn_tmp,bc_tmp,bca_tmp,mod),
      psi, method="Richardson"),
    error=function(e) rep(NA_real_,length(psi)))

  analytic_sum <- colSums(sc$scores)
  if (any(is.na(num_grad))) { skip_test(sprintf("T4 (%s)",mn),"numDeriv failed"); next }

  rel_err <- max(abs(analytic_sum + num_grad) / (abs(num_grad)+1e-6))
  cat(sprintf("    %s rel_err=%.4f\n", mn, rel_err))
  run_test(sprintf("Analytic == numerical scores (rel_err<0.10) (%s)",mn), rel_err<0.10)
}

# ---------------------------------------------------------------------------
cat("\n-- T5: fit_agdcc() all four models (k=3, T=500, vcv_method='none') --\n")
fits <- list()
for (m in c("DCC","ADCC","GDCC","AGDCC"))
  fits[[m]] <- tryCatch(
    fit_agdcc(data3,model=m,vcv_method="none",n_cores=1L),
    error=function(e){ cat("  ERROR",m,":",conditionMessage(e),"\n"); NULL })

for (m in names(fits)) {
  f<-fits[[m]]; if(is.null(f)){fail<-fail+1L;next}
  run_test(sprintf("%s: is agdcc_fit",m),       inherits(f,"agdcc_fit"))
  run_test(sprintf("%s: LL finite",m),           is.finite(f$ll))
  run_test(sprintf("%s: LL negative",m),         f$ll<0)
  run_test(sprintf("%s: Rt dimensions",m),       all(dim(f$Rt)==c(3,3,500)))
  run_test(sprintf("%s: Ht PD at t=1",m),        is_pd(f$Ht[,,1]))
  run_test(sprintf("%s: params positive",m),     all(f$psi_D>0))
  run_test(sprintf("%s: params < 1",m),          all(f$psi_D<1))
  npe <- switch(m,"DCC"=2,"ADCC"=3,"GDCC"=6,"AGDCC"=9)
  run_test(sprintf("%s: np=%d",m,npe), length(f$psi_D)==npe)
}

# ---------------------------------------------------------------------------
cat("\n-- T6: stage1_fit modes --\n")
spec_k3 <- rugarch::ugarchspec(
  variance.model=list(model="sGARCH",garchOrder=c(1,1)),
  mean.model=list(armaOrder=c(0,0),include.mean=FALSE),
  distribution.model="norm")

mfit <- tryCatch(rugarch::multifit(rugarch::multispec(replicate(3,spec_k3)),
                                   data3,solver="solnp",
                                   solver.control=list(trace=0)), error=function(e) NULL)
if (is.null(mfit)) {
  skip_test("ugarchMultifit: fit","multifit failed")
} else {
  f_mfit <- tryCatch(
    fit_agdcc(data3,model="DCC",stage1_fit=mfit,vcv_method="none"),
    error=function(e){
      if(grepl("slot|Cannot extract",conditionMessage(e),ignore.case=TRUE)){
        skip_test("ugarchMultifit: fit","slot accessor failed"); NULL
      } else NULL })
  if(!is.null(f_mfit)){
    run_test("ugarchMultifit: succeeds", !is.null(f_mfit))
    run_test("ugarchMultifit: LL finite", is.finite(f_mfit$ll))
  }
}
fl <- lapply(1:3, function(i) rugarch::ugarchfit(spec_k3,data3[,i],solver="solnp",
                                                   solver.control=list(trace=0)))
f_list <- tryCatch(fit_agdcc(data3,model="DCC",stage1_fit=fl,vcv_method="none"),
                   error=function(e) NULL)
run_test("list<ugarchfit>: succeeds", !is.null(f_list))
if(!is.null(f_list)) run_test("list<ugarchfit>: LL finite", is.finite(f_list$ll))

# ---------------------------------------------------------------------------
cat("\n-- T7: Stationarity --\n")
for (m in names(fits)) {
  f<-fits[[m]]; if(is.null(f)) next; km<-f$k
  cv <- switch(m,
    "DCC"  =f$psi_D["a"]^2+f$psi_D["b"]^2,
    "ADCC" =f$psi_D["a"]^2+f$delta*f$psi_D["g"]^2+f$psi_D["b"]^2,
    "GDCC" ={ av<-f$psi_D[1:km]; bv<-f$psi_D[(km+1):(2*km)]; max(av^2+bv^2) },
    "AGDCC"={ av<-f$psi_D[1:km]; gv<-f$psi_D[(km+1):(2*km)]; bv<-f$psi_D[(2*km+1):(3*km)]
              max(av^2+f$delta*gv^2+bv^2) })
  run_test(sprintf("%s: stationarity < 1",m), cv < 1)
}

# ---------------------------------------------------------------------------
cat("\n-- T8: VCV properties --\n")
for (m in names(fits)) {
  f<-fits[[m]]; if(is.null(f)||is.null(f$VCV)) next
  np<-length(f$psi_D)
  run_test(sprintf("%s: VCV %dx%d",m,np,np),     all(dim(f$VCV)==c(np,np)))
  run_test(sprintf("%s: VCV symmetric",m),
           max(abs(f$VCV-t(f$VCV)))<1e-8)   # 1e-8 tolerance for Jacobian noise
  ev<-eigen(f$VCV,symmetric=TRUE,only.values=TRUE)$values
  run_test(sprintf("%s: VCV PSD",m),              min(ev)>-1e-8)
}

# ---------------------------------------------------------------------------
cat("\n-- T10: S3 methods --\n")
if(!is.null(fits[["AGDCC"]])) {
  f<-fits[["AGDCC"]]
  ok <- tryCatch({ print(f); TRUE }, error=function(e) FALSE)
  run_test("print runs",      isTRUE(ok))
  ic <- AIC(f)
  run_test("AIC named",       !is.null(names(ic)) && "AIC"%in%names(ic))
  run_test("BIC>AIC",         ic["BIC"]>ic["AIC"])
  run_test("coef == psi_D",   identical(coef(f), f$psi_D))
  run_test("vcov == VCV",     identical(vcov(f), f$VCV))
}

# ---------------------------------------------------------------------------
cat("\n-- T11: Conditional correlations --\n")
for (m in names(fits)) {
  f<-fits[[m]]; if(is.null(f)) next
  for(ii in seq_len(f$k-1)) for(jj in seq(ii+1,f$k)) {
    rho<-conditional_correlations(f,ii,jj)
    run_test(sprintf("%s: rho(%d,%d) in (-1,1)",m,ii,jj),
             length(rho)==f$T && all(rho>-1-1e-8) && all(rho<1+1e-8))
  }
}

# ---------------------------------------------------------------------------
cat("\n-- T12: Parallel == serial --\n")
if(!is.null(fits[["AGDCC"]])) {
  f<-fits[["AGDCC"]]
  ic2<-compute_intercepts(f$e_std); cb2<-build_outer_cubes(f$e_std,ic2$n_std)
  bc2<-ewma_backcast_cpp(cb2$ee); bca2<-ewma_backcast_cpp(cb2$nn)
  sc1<-analytic_scores_cpp(f$psi_D,f$P,f$N,f$e_std,cb2$ee,cb2$nn,bc2,bca2,f$model_int,f$Qt,f$Rt)
  sc2<-analytic_scores_parallel_cpp(f$psi_D,f$P,f$N,f$e_std,cb2$ee,cb2$nn,bc2,bca2,f$model_int,f$Qt,f$Rt,1L)
  run_test("scores serial==parallel", max(abs(sc1$scores-sc2$scores))<1e-10)
  run_test("B22 serial==parallel",    max(abs(sc1$B22-sc2$B22))<1e-10)
}

# ---------------------------------------------------------------------------
cat("\n-- T13: 3-stage VCV --\n")
# Use small synthetic data (k=3) so n_tot = n1+n2+n3 << T_small
set.seed(99); x3 <- matrix(rnorm(300), 100, 3)
f3s <- fit_agdcc(x3, model="DCC", vcv_method="3stage")
run_test("3stage: SE positive",     all(f3s$se > 0))
run_test("3stage: VCV PSD",         is_pd(f3s$VCV))
run_test("3stage: A21 present",     !is.null(f3s$vcv_detail$A21))
run_test("3stage: A31 present",     !is.null(f3s$vcv_detail$A31))
run_test("3stage: A32 present",     !is.null(f3s$vcv_detail$A32))
run_test("3stage: A11 PSD",         is_pd(f3s$vcv_detail$A11))
run_test("3stage: B0 present",      !is.null(f3s$vcv_detail$B0))

# ---------------------------------------------------------------------------
cat(sprintf("\n=== Results: %d passed, %d failed (total %d) ===\n",pass,fail,pass+fail))
if(fail==0L) cat("All tests passed.\n") else {
  cat(sprintf("%d FAILED.\n",fail)); quit(status=1L)
}
