// agdcc_core.cpp
//
// C++ kernels for the AG-DCC estimator. Compiled via Rcpp::sourceCpp.
// Model codes: 1=DCC  2=ADCC  3=G-DCC  4=AG-DCC

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static inline mat normalise_qt(const mat& Q, vec& qsi) {
  const int k = (int)Q.n_rows;
  qsi.set_size(k);
  for (int i = 0; i < k; ++i)
    qsi(i) = 1.0 / std::sqrt(std::max(Q(i,i), 1e-14));
  mat R = Q % (qsi * qsi.t());
  R.diag().ones();
  return R;
}

static inline mat nearest_pd(const mat& M) {
  mat S = 0.5*(M + M.t());
  vec ev; mat ev_vec;
  eig_sym(ev, ev_vec, S);
  ev.transform([](double v){ return std::max(v, 1e-8); });
  return ev_vec * diagmat(ev) * ev_vec.t();
}

static inline bool safe_chol(const mat& R, mat& L) {
  if (chol(L, R, "lower")) return true;
  return chol(L, nearest_pd(R), "lower");
}

static inline double logdet_chol(const mat& L) {
  return 2.0 * accu(log(L.diag()));
}

static inline vec row_to_vec(const mat& M, int t) {
  return conv_to<vec>::from(M.row(t));
}

// Stationarity: spectral radius of (A⊗A + B⊗B) applied to Qt map.
// Scalar: a² + b² < 1.
// Diagonal: max_i(a_i² + b_i²) < 1  [NOT sum — sum grows with k].
// delta=1 is conservative; callers pass exact delta.
static inline double stat_val(const vec& p, int model, int k, double delta) {
  switch (model) {
    case 1: return p(0)*p(0) + p(1)*p(1);
    case 2: return p(0)*p(0) + delta*p(1)*p(1) + p(2)*p(2);
    case 3: { double mx=0; for(int i=0;i<k;++i){ double v=p(i)*p(i)+p(k+i)*p(k+i); if(v>mx) mx=v; } return mx; }
    default:{ double mx=0; for(int i=0;i<k;++i){ double v=p(i)*p(i)+delta*p(k+i)*p(k+i)+p(2*k+i)*p(2*k+i); if(v>mx) mx=v; } return mx; }
  }
}

// =============================================================================
// 1.  Qt / Rt RECURSION
// =============================================================================

// [[Rcpp::export]]
List qt_recursion_cpp(const arma::vec& psi_D,
                      const arma::mat& P,
                      const arma::mat& N,
                      const arma::cube& ee_arr,
                      const arma::cube& nn_arr,
                      const arma::mat& backcast,
                      const arma::mat& backcastAsym,
                      int model) {
  const int k=(int)P.n_rows, T=(int)ee_arr.n_slices;
  vec av(k), gv(k,fill::zeros), bv(k);
  switch(model){
    case 1: av.fill(psi_D(0)); bv.fill(psi_D(1)); break;
    case 2: av.fill(psi_D(0)); gv.fill(psi_D(1)); bv.fill(psi_D(2)); break;
    case 3: av=psi_D.subvec(0,k-1); bv=psi_D.subvec(k,2*k-1); break;
    default:av=psi_D.subvec(0,k-1);gv=psi_D.subvec(k,2*k-1);bv=psi_D.subvec(2*k,3*k-1);
  }
  mat AA=av*av.t(), GG=gv*gv.t(), BB=bv*bv.t();
  mat C = P - (AA%P) - (BB%P);
  if (any(gv>0.0)) C -= GG%N;

  cube Qt_out(k,k,T), Rt_out(k,k,T);
  mat Qt_prev=backcast;
  for(int t=0;t<T;++t){
    const mat& ee_lag=(t>0)?ee_arr.slice(t-1):backcast;
    const mat& nn_lag=(t>0)?nn_arr.slice(t-1):backcastAsym;
    mat Qt=C+(AA%ee_lag)+(GG%nn_lag)+(BB%Qt_prev);
    Qt=0.5*(Qt+Qt.t());
    Qt_out.slice(t)=Qt; Qt_prev=Qt;
    vec qi; Rt_out.slice(t)=normalise_qt(Qt,qi);
  }
  return List::create(Named("Qt")=Qt_out, Named("Rt")=Rt_out);
}

// =============================================================================
// 2.  STAGE-3 LOG-LIKELIHOOD
//
// stage3_qll_fast_cpp: single-pass recursion + LL without storing Qt/Rt cubes.
//   Avoids allocating 2*k*k*T doubles on every objective call. Used by
//   stage3_obj and stage3_qll_grid_cpp to keep Stage-3 optimization fast even
//   for large k (k=34 saves ~25 MB per call and eliminates GC pressure).
//
// stage3_qll_cpp: legacy version (stores full Qt/Rt cubes) retained for the
//   qt_recursion → analytic_scores two-step path used in post-estimation VCV.
// =============================================================================

// [[Rcpp::export]]
double stage3_qll_fast_cpp(const arma::vec& psi_D,
                            const arma::mat& P,
                            const arma::mat& N,
                            const arma::mat& e_std,
                            const arma::cube& ee_arr,
                            const arma::cube& nn_arr,
                            const arma::mat& backcast,
                            const arma::mat& backcastAsym,
                            int model) {
  const int k=(int)P.n_rows, T=(int)e_std.n_rows;
  for(uword i=0;i<psi_D.n_elem;++i)
    if(psi_D(i)<=0.0||psi_D(i)>=1.0) return 1e9;
  if(stat_val(psi_D,model,k,1.0)>=0.9999) return 1e9;

  vec av(k),gv(k,fill::zeros),bv(k);
  switch(model){
    case 1: av.fill(psi_D(0)); bv.fill(psi_D(1)); break;
    case 2: av.fill(psi_D(0)); gv.fill(psi_D(1)); bv.fill(psi_D(2)); break;
    case 3: av=psi_D.subvec(0,k-1); bv=psi_D.subvec(k,2*k-1); break;
    default:av=psi_D.subvec(0,k-1);gv=psi_D.subvec(k,2*k-1);bv=psi_D.subvec(2*k,3*k-1);
  }
  mat AA=av*av.t(), GG=gv*gv.t(), BB=bv*bv.t();
  mat C=P-(AA%P)-(BB%P);
  if(any(gv>0.0)) C-=GG%N;

  double ll=0.0;
  mat Qt_prev=backcast;
  for(int t=0;t<T;++t){
    const mat& eel=(t>0)?ee_arr.slice(t-1):backcast;
    const mat& nnl=(t>0)?nn_arr.slice(t-1):backcastAsym;
    mat Qt=C+(AA%eel)+(GG%nnl)+(BB%Qt_prev);
    Qt=0.5*(Qt+Qt.t());
    vec qi; mat Rt=normalise_qt(Qt,qi);
    mat L;
    if(!safe_chol(Rt,L)){ll+=1e6; Qt_prev=Qt; continue;}
    vec et=row_to_vec(e_std,t);
    vec z=solve(trimatl(L),et);
    ll+=0.5*(logdet_chol(L)+dot(z,z));
    Qt_prev=Qt;
  }
  return std::isfinite(ll)?ll:1e9;
}
//
// Stationarity gate: >= 0.9999  (SOFT, not hard 1.0)
// Reason: numDeriv Richardson extrapolation steps h ~ 1e-4 * |psi|.
// At psi=(0.05, 0.90), stat_val = 0.8125. A step of +1e-4 on b gives
// stat_val(0.90001) = 0.8101 < 0.9999 — safe. Using hard gate 1.0 instead
// would require stat_val to exceed 1.0, which only happens if a parameter
// hits near 1.0 individually, still caught by the box check psi(i) >= 1.0.
// The soft gate catches the penalty region consistently with the R-side penalty.

// [[Rcpp::export]]
double stage3_qll_cpp(const arma::vec& psi_D,
                      const arma::mat& P,
                      const arma::mat& N,
                      const arma::mat& e_std,
                      const arma::cube& ee_arr,
                      const arma::cube& nn_arr,
                      const arma::mat& backcast,
                      const arma::mat& backcastAsym,
                      int model) {
  const int k=(int)P.n_rows, T=(int)e_std.n_rows;
  for(uword i=0;i<psi_D.n_elem;++i)
    if(psi_D(i)<=0.0||psi_D(i)>=1.0) return 1e9;
  // Soft stationarity gate (delta=1 conservative)
  if(stat_val(psi_D,model,k,1.0)>=0.9999) return 1e9;

  List rec;
  try{ rec=qt_recursion_cpp(psi_D,P,N,ee_arr,nn_arr,backcast,backcastAsym,model); }
  catch(...){ return 1e9; }
  const cube& Rt_arr=rec["Rt"];

  double ll=0.0;
  for(int t=0;t<T;++t){
    mat Rt=Rt_arr.slice(t); mat L;
    if(!safe_chol(Rt,L)){ll+=1e6;continue;}
    vec et=row_to_vec(e_std,t);
    vec z=solve(trimatl(L),et);
    ll+=0.5*(logdet_chol(L)+dot(z,z));
  }
  return std::isfinite(ll)?ll:1e9;
}

// =============================================================================
// 3.  ANALYTIC DING-ENGLE SCORE RECURSION  (serial)
// =============================================================================

// [[Rcpp::export]]
List analytic_scores_cpp(const arma::vec& psi_D,
                         const arma::mat& P,
                         const arma::mat& N,
                         const arma::mat& e_std,
                         const arma::cube& ee_arr,
                         const arma::cube& nn_arr,
                         const arma::mat& backcast,
                         const arma::mat& backcastAsym,
                         int model,
                         const arma::cube& Qt_arr,
                         const arma::cube& Rt_arr) {
  const int k=(int)P.n_rows,T=(int)e_std.n_rows,np=(int)psi_D.n_elem;
  vec av(k),gv(k,fill::zeros),bv(k);
  switch(model){
    case 1: av.fill(psi_D(0)); bv.fill(psi_D(1)); break;
    case 2: av.fill(psi_D(0)); gv.fill(psi_D(1)); bv.fill(psi_D(2)); break;
    case 3: av=psi_D.subvec(0,k-1); bv=psi_D.subvec(k,2*k-1); break;
    default:av=psi_D.subvec(0,k-1);gv=psi_D.subvec(k,2*k-1);bv=psi_D.subvec(2*k,3*k-1);
  }
  mat BB=bv*bv.t();
  mat sc(T,np,fill::zeros);
  std::vector<mat> dQp(np,zeros<mat>(k,k)), dQc(np,zeros<mat>(k,k));

  for(int t=0;t<T;++t){
    const mat& Qt=Qt_arr.slice(t);
    const mat& Rt=Rt_arr.slice(t);
    vec qi(k);
    for(int i=0;i<k;++i) qi(i)=1.0/std::sqrt(std::max(Qt(i,i),1e-14));
    mat Qsi=diagmat(qi);
    mat Ri; if(!inv_sympd(Ri,Rt)) Ri=pinv(Rt);
    vec et=row_to_vec(e_std,t), ut=Ri*et;
    mat Om=Ri-ut*ut.t();
    const mat& eel=(t>0)?ee_arr.slice(t-1):backcast;
    const mat& nnl=(t>0)?nn_arr.slice(t-1):backcastAsym;
    const mat& Ql=(t>0)?Qt_arr.slice(t-1):backcast;

    // dQt/dpsi: correct chain-rule derivatives of Qt = P(I-aa'-bb') - N*gg' + aa'@ee + gg'@nn + bb'@Qt-1
    // For scalar (model 1,2): d(aa')/da = 2a, d(bb')/db = 2b, persistence uses b^2
    // For diagonal (model 3,4): d(aa')/da_i = e_i*a' + a*e_i' (NOT 2*a_i*S_i)
    double b2sq=(model<=2)?psi_D(model==1?1:2)*psi_D(model==1?1:2):0.0;
    switch(model){
      case 1:
        dQc[0]=2.0*psi_D(0)*(eel-P)+b2sq*dQp[0];
        dQc[1]=2.0*psi_D(1)*(Ql-P) +b2sq*dQp[1]; break;
      case 2:
        dQc[0]=2.0*psi_D(0)*(eel-P)+b2sq*dQp[0];
        dQc[1]=2.0*psi_D(1)*(nnl-N)+b2sq*dQp[1];
        dQc[2]=2.0*psi_D(2)*(Ql-P) +b2sq*dQp[2]; break;
      case 3:
        for(int i=0;i<k;++i){
          mat dA(k,k,fill::zeros); dA.row(i)=av.t(); dA.col(i)+=av;
          dQc[i]  =dA%(eel-P)+BB%dQp[i];
          mat dB(k,k,fill::zeros); dB.row(i)=bv.t(); dB.col(i)+=bv;
          dQc[k+i]=dB%(Ql-P)+BB%dQp[k+i];
        } break;
      default:
        for(int i=0;i<k;++i){
          mat dA(k,k,fill::zeros); dA.row(i)=av.t(); dA.col(i)+=av;
          dQc[i]    =dA%(eel-P)+BB%dQp[i];
          mat dG(k,k,fill::zeros); dG.row(i)=gv.t(); dG.col(i)+=gv;
          dQc[k+i]  =dG%(nnl-N)+BB%dQp[k+i];
          mat dB(k,k,fill::zeros); dB.row(i)=bv.t(); dB.col(i)+=bv;
          dQc[2*k+i]=dB%(Ql-P)+BB%dQp[2*k+i];
        }
    }
    for(int j=0;j<np;++j){
      const mat& dQ=dQc[j];
      mat ddg=diagmat(dQ.diag());
      mat dR=Qsi*dQ*Qsi - 0.5*(Qsi*ddg*Qsi*Rt + Rt*Qsi*ddg*Qsi);
      sc(t,j)=-0.5*accu(Om%dR);
    }
    dQp=dQc;
  }
  mat B22=sc.t()*sc/(double)T;
  return List::create(Named("scores")=sc, Named("B22")=B22);
}

// =============================================================================
// 4.  PARALLEL SCORE RECURSION
// =============================================================================

// [[Rcpp::export]]
List analytic_scores_parallel_cpp(const arma::vec& psi_D,
                                   const arma::mat& P,
                                   const arma::mat& N,
                                   const arma::mat& e_std,
                                   const arma::cube& ee_arr,
                                   const arma::cube& nn_arr,
                                   const arma::mat& backcast,
                                   const arma::mat& backcastAsym,
                                   int model,
                                   const arma::cube& Qt_arr,
                                   const arma::cube& Rt_arr,
                                   int n_cores=1){
#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#endif
  const int k=(int)P.n_rows,T=(int)e_std.n_rows,np=(int)psi_D.n_elem;
  vec av(k),gv(k,fill::zeros),bv(k);
  switch(model){
    case 1: av.fill(psi_D(0)); bv.fill(psi_D(1)); break;
    case 2: av.fill(psi_D(0)); gv.fill(psi_D(1)); bv.fill(psi_D(2)); break;
    case 3: av=psi_D.subvec(0,k-1); bv=psi_D.subvec(k,2*k-1); break;
    default:av=psi_D.subvec(0,k-1);gv=psi_D.subvec(k,2*k-1);bv=psi_D.subvec(2*k,3*k-1);
  }
  mat BB=bv*bv.t();
  std::vector<double> sr((size_t)T*np,0.0);
  std::vector<mat> dQp(np,zeros<mat>(k,k)),dQc(np,zeros<mat>(k,k));

  for(int t=0;t<T;++t){
    const mat& Qt=Qt_arr.slice(t);
    const mat& Rt=Rt_arr.slice(t);
    vec qi(k);
    for(int i=0;i<k;++i) qi(i)=1.0/std::sqrt(std::max(Qt(i,i),1e-14));
    mat Qsi=diagmat(qi);
    mat Ri; if(!inv_sympd(Ri,Rt)) Ri=pinv(Rt);
    vec et=row_to_vec(e_std,t), ut=Ri*et;
    mat Om=Ri-ut*ut.t();
    const mat& eel=(t>0)?ee_arr.slice(t-1):backcast;
    const mat& nnl=(t>0)?nn_arr.slice(t-1):backcastAsym;
    const mat& Ql=(t>0)?Qt_arr.slice(t-1):backcast;

    double b2sq_p=(model<=2)?psi_D(model==1?1:2)*psi_D(model==1?1:2):0.0;
    switch(model){
      case 1:
        dQc[0]=2.0*psi_D(0)*(eel-P)+b2sq_p*dQp[0];
        dQc[1]=2.0*psi_D(1)*(Ql-P) +b2sq_p*dQp[1]; break;
      case 2:
        dQc[0]=2.0*psi_D(0)*(eel-P)+b2sq_p*dQp[0];
        dQc[1]=2.0*psi_D(1)*(nnl-N)+b2sq_p*dQp[1];
        dQc[2]=2.0*psi_D(2)*(Ql-P) +b2sq_p*dQp[2]; break;
      case 3:
        for(int i=0;i<k;++i){
          mat dA(k,k,fill::zeros); dA.row(i)=av.t(); dA.col(i)+=av;
          dQc[i]  =dA%(eel-P)+BB%dQp[i];
          mat dB(k,k,fill::zeros); dB.row(i)=bv.t(); dB.col(i)+=bv;
          dQc[k+i]=dB%(Ql-P)+BB%dQp[k+i];
        } break;
      default:
        for(int i=0;i<k;++i){
          mat dA(k,k,fill::zeros); dA.row(i)=av.t(); dA.col(i)+=av;
          dQc[i]    =dA%(eel-P)+BB%dQp[i];
          mat dG(k,k,fill::zeros); dG.row(i)=gv.t(); dG.col(i)+=gv;
          dQc[k+i]  =dG%(nnl-N)+BB%dQp[k+i];
          mat dB(k,k,fill::zeros); dB.row(i)=bv.t(); dB.col(i)+=bv;
          dQc[2*k+i]=dB%(Ql-P)+BB%dQp[2*k+i];
        }
    }
    int npl=np; double* srp=sr.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int j=0;j<npl;++j){
      const mat& dQ=dQc[j];
      mat ddg=diagmat(dQ.diag());
      mat dR=Qsi*dQ*Qsi-0.5*(Qsi*ddg*Qsi*Rt+Rt*Qsi*ddg*Qsi);
      srp[(size_t)t+(size_t)j*T]=-0.5*accu(Om%dR);
    }
    dQp=dQc;
  }
  mat sc(sr.data(),T,np);
  mat B22=sc.t()*sc/(double)T;
  return List::create(Named("scores")=sc, Named("B22")=B22);
}

// =============================================================================
// 5.  EWMA BACK-CAST
// =============================================================================

// ns = sqrt(T) matches the MFE toolbox convention for the EWMA burn-in window.
// [[Rcpp::export]]
arma::mat ewma_backcast_cpp(const arma::cube& X, double lambda=0.94){
  const int k=(int)X.n_rows, ns=(int)std::sqrt((double)X.n_slices);
  mat BC(k,k,fill::zeros);
  double cum=0.0; std::vector<double> ww(ns);
  for(int j=0;j<ns;++j){ ww[j]=(1.0-lambda)*std::pow(lambda,j); cum+=ww[j]; }
  for(int j=0;j<ns;++j) BC+=(ww[j]/cum)*X.slice(j);
  return BC;
}

// =============================================================================
// 7.  GRID EVALUATION
// =============================================================================

// [[Rcpp::export]]
arma::vec stage3_qll_grid_cpp(const arma::mat& grid,
                               const arma::mat& P,
                               const arma::mat& N,
                               const arma::mat& e_std,
                               const arma::cube& ee_arr,
                               const arma::cube& nn_arr,
                               const arma::mat& backcast,
                               const arma::mat& backcastAsym,
                               int model, int n_cores=1){
#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#endif
  const int ng=(int)grid.n_rows; vec ll(ng);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for(int g=0;g<ng;++g){
    vec p=conv_to<vec>::from(grid.row(g));
    ll(g)=stage3_qll_fast_cpp(p,P,N,e_std,ee_arr,nn_arr,backcast,backcastAsym,model);
  }
  return ll;
}

// =============================================================================
// 8.  BUILD OUTER CUBES  (replaces R loop — avoids GC pressure)
// =============================================================================

// [[Rcpp::export]]
List build_outer_cubes_cpp(const arma::mat& e_std, const arma::mat& n_std) {
  const int T = (int)e_std.n_rows, k = (int)e_std.n_cols;
  cube ee(k, k, T), nn(k, k, T);
  for (int t = 0; t < T; ++t) {
    vec et = conv_to<vec>::from(e_std.row(t));
    vec nt = conv_to<vec>::from(n_std.row(t));
    ee.slice(t) = et * et.t();
    nn.slice(t) = nt * nt.t();
  }
  return List::create(Named("ee") = ee, Named("nn") = nn);
}

// =============================================================================
// 10. SINGLE-THREADED Qt RECURSION + MEAN SCORES (internal helper)
//
// Identical logic to qt_and_mean_scores_cpp but without omp_set_num_threads
// and without the inner OpenMP pragma. Safe to call from OMP worker threads
// in compute_A32_cpp (each thread needs its own independent Qt recursion).
// Requires a thread-safe BLAS (OpenBLAS/MKL, standard on Windows R builds).
// =============================================================================

static arma::vec qt_mean_scores_1t(const arma::vec& psi_D,
                                    const arma::mat& P, const arma::mat& N,
                                    const arma::mat& e_std,
                                    const arma::cube& ee, const arma::cube& nn,
                                    const arma::mat& bc, const arma::mat& bca,
                                    int model) {
  const int k = (int)P.n_rows, T = (int)e_std.n_rows, np = (int)psi_D.n_elem;
  vec av(k), gv(k, fill::zeros), bv(k);
  switch (model) {
    case 1: av.fill(psi_D(0)); bv.fill(psi_D(1)); break;
    case 2: av.fill(psi_D(0)); gv.fill(psi_D(1)); bv.fill(psi_D(2)); break;
    case 3: av = psi_D.subvec(0,k-1); bv = psi_D.subvec(k,2*k-1); break;
    default: av = psi_D.subvec(0,k-1); gv = psi_D.subvec(k,2*k-1); bv = psi_D.subvec(2*k,3*k-1);
  }
  mat AA = av*av.t(), GG = gv*gv.t(), BB = bv*bv.t();
  mat C = P - (AA%P) - (BB%P);
  if (any(gv > 0.0)) C -= GG%N;

  vec mean_sc(np, fill::zeros);
  std::vector<mat> dQp(np, zeros<mat>(k,k)), dQc(np, zeros<mat>(k,k));
  mat Qt_prev = bc;

  for (int t = 0; t < T; ++t) {
    const mat& eel = (t > 0) ? ee.slice(t-1) : bc;
    const mat& nnl = (t > 0) ? nn.slice(t-1) : bca;
    const mat& Ql  = Qt_prev;
    mat Qt = C + (AA%eel) + (GG%nnl) + (BB%Ql);
    Qt = 0.5*(Qt + Qt.t());
    vec qi; mat Rt = normalise_qt(Qt, qi);
    mat Qsi = diagmat(qi);
    mat Ri; if (!inv_sympd(Ri, Rt)) Ri = pinv(Rt);
    vec et = conv_to<vec>::from(e_std.row(t));
    mat Om = Ri - (Ri*et)*(Ri*et).t();
    double b2sq = (model<=2) ? psi_D(model==1?1:2)*psi_D(model==1?1:2) : 0.0;
    switch (model) {
      case 1:
        dQc[0] = 2.0*psi_D(0)*(eel-P) + b2sq*dQp[0];
        dQc[1] = 2.0*psi_D(1)*(Ql -P) + b2sq*dQp[1]; break;
      case 2:
        dQc[0] = 2.0*psi_D(0)*(eel-P) + b2sq*dQp[0];
        dQc[1] = 2.0*psi_D(1)*(nnl-N) + b2sq*dQp[1];
        dQc[2] = 2.0*psi_D(2)*(Ql -P) + b2sq*dQp[2]; break;
      case 3:
        for (int i = 0; i < k; ++i) {
          mat dA(k,k,fill::zeros); dA.row(i)=av.t(); dA.col(i)+=av;
          dQc[i]   = dA%(eel-P) + BB%dQp[i];
          mat dB(k,k,fill::zeros); dB.row(i)=bv.t(); dB.col(i)+=bv;
          dQc[k+i] = dB%(Ql-P)  + BB%dQp[k+i];
        } break;
      default:
        for (int i = 0; i < k; ++i) {
          mat dA(k,k,fill::zeros); dA.row(i)=av.t(); dA.col(i)+=av;
          dQc[i]     = dA%(eel-P) + BB%dQp[i];
          mat dG(k,k,fill::zeros); dG.row(i)=gv.t(); dG.col(i)+=gv;
          dQc[k+i]   = dG%(nnl-N) + BB%dQp[k+i];
          mat dB(k,k,fill::zeros); dB.row(i)=bv.t(); dB.col(i)+=bv;
          dQc[2*k+i] = dB%(Ql-P)  + BB%dQp[2*k+i];
        }
    }
    for (int j = 0; j < np; ++j) {
      const mat& dQ = dQc[j];
      mat ddg = diagmat(dQ.diag());
      mat dR = Qsi*dQ*Qsi - 0.5*(Qsi*ddg*Qsi*Rt + Rt*Qsi*ddg*Qsi);
      mean_sc(j) += -0.5*accu(Om%dR);
    }
    Qt_prev = Qt; dQp = dQc;
  }
  mean_sc /= (double)T;
  return mean_sc;
}

// =============================================================================
// 11. COMPUTE A32 IN C++ WITH OPENMP OVER PAIRS
//
// Replaces the R-level compute_A32 (parLapply over pairs). Benefits:
//   - No cluster setup/teardown overhead
//   - Cubes (ee, nn, bc, bca) allocated once and shared read-only across threads
//   - Each thread: local copy of P_p or N_p (stack), independent Qt recursion
//   - schedule(dynamic): load-balances R-bar vs N-bar pairs
//
// Pair ordering (column-major, 0-based): R-bar block j=0..k-2, i=j+1..k-1;
// N-bar block (ADCC/AGDCC) j=0..k-1, i=j..k-1. Matches .pairs_lower / .pairs_lower_diag.
// =============================================================================

// [[Rcpp::export]]
arma::mat compute_A32_cpp(const arma::vec& psi_D,
                           const arma::mat& P, const arma::mat& N,
                           const arma::mat& e_std,
                           const arma::cube& ee, const arma::cube& nn,
                           const arma::mat& bc, const arma::mat& bca,
                           int model, int n_cores = 1, double h = 1e-5) {
  const int k = (int)P.n_rows, np = (int)psi_D.n_elem;
  const bool is_asym = (model == 2 || model == 4);
  const int n2_sym  = k*(k-1)/2;
  const int n2_asym = is_asym ? k*(k+1)/2 : 0;
  const int n2      = n2_sym + n2_asym;

  std::vector<int> ri(n2), ci(n2);
  int idx = 0;
  for (int j = 0; j < k-1; ++j)
    for (int i = j+1; i < k; ++i) { ri[idx]=i; ci[idx]=j; ++idx; }
  if (is_asym)
    for (int j = 0; j < k; ++j)
      for (int i = j; i < k; ++i) { ri[idx]=i; ci[idx]=j; ++idx; }

  vec s3_base = qt_mean_scores_1t(psi_D, P, N, e_std, ee, nn, bc, bca, model);
  mat A32(np, n2, fill::zeros);

#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#pragma omp parallel for schedule(dynamic)
#endif
  for (int p = 0; p < n2; ++p) {
    vec s3_p;
    if (p < n2_sym) {
      mat P_p = P;
      P_p(ri[p], ci[p]) += h; P_p(ci[p], ri[p]) += h;
      s3_p = qt_mean_scores_1t(psi_D, P_p, N, e_std, ee, nn, bc, bca, model);
    } else {
      mat N_p = N;
      N_p(ri[p], ci[p]) += h;
      if (ri[p] != ci[p]) N_p(ci[p], ri[p]) += h;
      s3_p = qt_mean_scores_1t(psi_D, P, N_p, e_std, ee, nn, bc, bca, model);
    }
    A32.col(p) = -(s3_p - s3_base) / h;
  }
  return A32;
}

// =============================================================================
// 9.  COMBINED Qt RECURSION + MEAN ANALYTIC SCORES
//     Returns only mean scores (length-np vec) without passing Qt/Rt cubes
//     back through R. Used inside numDeriv::jacobian in vcv_mfe to avoid
//     repeated large-cube allocation/GC across ~4*np evaluations.
// =============================================================================

// [[Rcpp::export]]
arma::vec qt_and_mean_scores_cpp(const arma::vec& psi_D,
                                  const arma::mat& P,
                                  const arma::mat& N,
                                  const arma::mat& e_std,
                                  const arma::cube& ee_arr,
                                  const arma::cube& nn_arr,
                                  const arma::mat& backcast,
                                  const arma::mat& backcastAsym,
                                  int model,
                                  int n_cores = 1) {
#ifdef _OPENMP
  omp_set_num_threads(n_cores);
#endif
  const int k = (int)P.n_rows, T = (int)e_std.n_rows, np = (int)psi_D.n_elem;

  vec av(k), gv(k, fill::zeros), bv(k);
  switch (model) {
    case 1: av.fill(psi_D(0)); bv.fill(psi_D(1)); break;
    case 2: av.fill(psi_D(0)); gv.fill(psi_D(1)); bv.fill(psi_D(2)); break;
    case 3: av = psi_D.subvec(0,k-1); bv = psi_D.subvec(k,2*k-1); break;
    default: av = psi_D.subvec(0,k-1); gv = psi_D.subvec(k,2*k-1); bv = psi_D.subvec(2*k,3*k-1);
  }
  mat AA = av*av.t(), GG = gv*gv.t(), BB = bv*bv.t();
  mat C  = P - (AA%P) - (BB%P);
  if (any(gv > 0.0)) C -= GG%N;

  vec mean_sc(np, fill::zeros);
  std::vector<mat> dQp(np, zeros<mat>(k,k)), dQc(np, zeros<mat>(k,k));
  mat Qt_prev = backcast;

  for (int t = 0; t < T; ++t) {
    const mat& eel = (t > 0) ? ee_arr.slice(t-1) : backcast;
    const mat& nnl = (t > 0) ? nn_arr.slice(t-1) : backcastAsym;
    const mat& Ql  = Qt_prev;

    mat Qt = C + (AA%eel) + (GG%nnl) + (BB%Ql);
    Qt = 0.5*(Qt + Qt.t());
    vec qi; mat Rt = normalise_qt(Qt, qi);
    mat Qsi = diagmat(qi);

    mat Ri; if (!inv_sympd(Ri, Rt)) Ri = pinv(Rt);
    vec et  = conv_to<vec>::from(e_std.row(t));
    mat Om  = Ri - (Ri*et)*(Ri*et).t();

    double b2sq = (model<=2) ? psi_D(model==1?1:2)*psi_D(model==1?1:2) : 0.0;
    switch (model) {
      case 1:
        dQc[0] = 2.0*psi_D(0)*(eel-P) + b2sq*dQp[0];
        dQc[1] = 2.0*psi_D(1)*(Ql -P) + b2sq*dQp[1]; break;
      case 2:
        dQc[0] = 2.0*psi_D(0)*(eel-P) + b2sq*dQp[0];
        dQc[1] = 2.0*psi_D(1)*(nnl-N) + b2sq*dQp[1];
        dQc[2] = 2.0*psi_D(2)*(Ql -P) + b2sq*dQp[2]; break;
      case 3:
        for (int i = 0; i < k; ++i) {
          mat dA(k,k,fill::zeros); dA.row(i)=av.t(); dA.col(i)+=av;
          dQc[i]   = dA%(eel-P) + BB%dQp[i];
          mat dB(k,k,fill::zeros); dB.row(i)=bv.t(); dB.col(i)+=bv;
          dQc[k+i] = dB%(Ql-P)  + BB%dQp[k+i];
        } break;
      default:
        for (int i = 0; i < k; ++i) {
          mat dA(k,k,fill::zeros); dA.row(i)=av.t(); dA.col(i)+=av;
          dQc[i]     = dA%(eel-P) + BB%dQp[i];
          mat dG(k,k,fill::zeros); dG.row(i)=gv.t(); dG.col(i)+=gv;
          dQc[k+i]   = dG%(nnl-N) + BB%dQp[k+i];
          mat dB(k,k,fill::zeros); dB.row(i)=bv.t(); dB.col(i)+=bv;
          dQc[2*k+i] = dB%(Ql-P)  + BB%dQp[2*k+i];
        }
    }

    int npl = np; double* ms = mean_sc.memptr();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int j = 0; j < npl; ++j) {
      const mat& dQ = dQc[j];
      mat ddg = diagmat(dQ.diag());
      mat dR  = Qsi*dQ*Qsi - 0.5*(Qsi*ddg*Qsi*Rt + Rt*Qsi*ddg*Qsi);
      ms[j]  += -0.5*accu(Om%dR);
    }

    Qt_prev = Qt;
    dQp = dQc;
  }
  mean_sc /= (double)T;
  return mean_sc;
}
