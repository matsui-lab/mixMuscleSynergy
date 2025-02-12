// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//----------------------------------------------------
// dmvnormFA_rowwiseCpp: (シングルスレッド) 多次元正規対数密度
//----------------------------------------------------

// [[Rcpp::export]]
arma::vec dmvnormFA_rowwiseCpp(const arma::mat &X,
                               const arma::vec &mu,
                               const arma::mat &Lambda,
                               const arma::vec &Psi_diag)
{
  int T = X.n_rows;
  int M = X.n_cols;

  // Sigma = Lambda * Lambda.t() + diag(Psi_diag)
  mat Sigma = Lambda * Lambda.t();
  for(int m=0; m<M; m++){
    Sigma(m,m) += Psi_diag[m];
  }
  // Numerical stabilization
  Sigma.diag() += 1e-8;

  // Chol(Sigma)
  mat L = chol(Sigma, "lower");
  double logdetSigma = 0.0;
  for(int m=0; m<M; m++){
    logdetSigma += 2.0 * std::log(L(m,m));
  }
  double log2pi = std::log(2.0 * M_PI);
  double Mlog2pi = M * log2pi;

  // L^-1
  mat Linv = inv(trimatl(L));

  // Output vector
  arma::vec out(T, fill::zeros);

  // Row-wise loop
  for(int i = 0; i < T; i++){
    vec diff = X.row(i).t() - mu;
    vec z = Linv * diff;
    double quadform = dot(z,z);
    double val = -0.5 * (Mlog2pi + logdetSigma + quadform);
    out[i] = val;
  }

  return out;
}

//----------------------------------------------------
// faEMupdateCpp: (単一クラスタ) FactorAnalyzer EM
//----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List faEMupdateCpp(const arma::mat &X, int r, int nIterFA)
{
  int n = X.n_rows;
  int M = X.n_cols;

  // Mean
  rowvec mu_r = mean(X, 0);
  mat Xc = X;
  for(int i=0; i<n; i++){
    Xc.row(i) -= mu_r;
  }

  // Initial Lambda => via PCA (SVD)
  mat U, V;
  vec S;
  svd(U, S, V, Xc);
  mat Vr = V.cols(0, r-1);
  vec sr = S.subvec(0, r-1);
  mat Lambda = Vr * diagmat(sr / std::sqrt((double)n));

  // Initial psi
  vec psi(M, fill::ones);
  double extraVar = 0.0;
  if(r < M){
    for(int j=r; j<M; j++){
      extraVar += S[j];
    }
    extraVar /= double(M - r) * double(n);
  } else {
    extraVar = 1.0;
  }
  psi *= extraVar;

  // EM iterations
  for(int it=0; it<nIterFA; it++){
    mat Sig = Lambda * Lambda.t();
    for(int m=0; m<M; m++){
      Sig(m,m) += psi[m];
    }
    // Numerical stabilization
    Sig.diag() += 1e-8;

    mat L = chol(Sig, "lower");
    mat invL = inv(trimatl(L));
    mat invSigma = invL.t() * invL;

    mat W = Lambda.t() * invSigma;
    mat A = eye(r, r) + W * Lambda;
    mat A_chol = chol(A, "lower");
    mat A_inv  = inv(trimatl(A_chol));
    mat Vz     = A_inv.t() * A_inv;

    mat WX = W * Xc.t();
    mat EZ = Vz * WX;

    mat EZZt = (EZ * EZ.t()) / double(n);
    for(int j=0; j<r; j++){
      EZZt(j,j) += Vz(j,j);
    }

    mat XcEZt = (Xc.t() * EZ.t()) / double(n);
    mat Lambda_new = XcEZt * inv(EZZt + 1e-12*eye(r,r));

    mat Sxx = (Xc.t() * Xc) / double(n);
    mat tmp = Lambda_new * XcEZt.t();
    mat residual = Sxx - tmp;
    for(int m=0; m<M; m++){
      double val = residual(m,m);
      if(val < 1e-12) val = 1e-12;
      psi[m] = val;
    }

    Lambda = Lambda_new;
  }

  // Return
  List out;
  out["mu"]     = mu_r.t();
  out["Lambda"] = Lambda;
  out["Psi"]    = diagmat(psi);
  return out;
}

//----------------------------------------------------
// mfaTimeseriesCpp: Mixture of Factor Analyzers (hard assignment)
//    - If z_init is provided, use it as the initial assignment
//    - Otherwise default to (i % K) + 1
//----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List mfaTimeseriesCpp(Rcpp::List list_of_data,
                            int K,
                            int r,
                            int max_iter,
                            int nIterFA,
                            double tol,
                            Rcpp::Nullable<Rcpp::IntegerVector> z_init = R_NilValue)
{
  int N = list_of_data.size();
  if(N < 1) {
    stop("No data in list_of_data");
  }

  // Get M from the first subject
  NumericMatrix firstMat = list_of_data[0];
  int M = firstMat.ncol();

  // Initialize cluster labels zR
  IntegerVector zR(N);

  // If user provided an initial assignment
  if(z_init.isNotNull()){
    Rcpp::IntegerVector zz(z_init.get());
    if(zz.size() != N){
      stop("z_init size != N");
    }
    // Check 1..K
    for(int i=0; i<N; i++){
      if(zz[i] < 1 || zz[i] > K){
        stop("z_init has out-of-range cluster label");
      }
      zR[i] = zz[i];
    }
  } else {
    // Default init
    for(int i=0; i<N; i++){
      zR[i] = (i % K) + 1;
    }
  }

  // Parameters
  std::vector<arma::mat> Lambda(K);
  std::vector<arma::vec> mu(K);
  std::vector<arma::mat> Psi(K);

  // Random init
  std::srand(123);
  for(int k=0; k<K; k++){
    Lambda[k] = 0.01 * arma::randn<arma::mat>(M, r);
    mu[k]     = arma::vec(M, fill::zeros);
    Psi[k]    = arma::eye<arma::mat>(M,M)*0.1;
  }

  arma::vec pi_k(K, fill::value(1.0 / double(K)));
  double old_loglik = -1e15;

  // EM loop
  for(int iterEM = 1; iterEM <= max_iter; iterEM++){
    //==============================
    // M-step
    //==============================
    for(int k=0; k<K; k++){
      // gather indices for cluster k
      std::vector<int> idx;
      idx.reserve(N);
      for(int i=0; i<N; i++){
        if(zR[i] == (k+1)){
          idx.push_back(i);
        }
      }
      if(idx.empty()) {
        // skip empty cluster
        continue;
      }

      // concatenate data
      long totalT = 0;
      for(size_t ii=0; ii<idx.size(); ii++){
        NumericMatrix mat_i = list_of_data[idx[ii]];
        totalT += mat_i.nrow();
      }
      arma::mat X_k(totalT, M, fill::zeros);
      {
        long rowpos = 0;
        for(size_t ii=0; ii<idx.size(); ii++){
          NumericMatrix mat_i = list_of_data[idx[ii]];
          int Ti = mat_i.nrow();
          for(int t=0; t<Ti; t++){
            for(int mm=0; mm<M; mm++){
              X_k(rowpos, mm) = mat_i(t, mm);
            }
            rowpos++;
          }
        }
      }

      // single-cluster Factor Analysis
      Rcpp::List fares = faEMupdateCpp(X_k, r, nIterFA);
      mu[k]     = as<arma::vec>(fares["mu"]);
      Lambda[k] = as<arma::mat>(fares["Lambda"]);
      Psi[k]    = as<arma::mat>(fares["Psi"]);
    }

    //==============================
    // E-step (hard assignment)
    //==============================
    arma::mat logLik_i(N, K, fill::zeros);
    arma::vec logPi(K);
    for(int k2=0; k2<K; k2++){
      logPi[k2] = std::log(pi_k[k2] + 1e-16);
    }

    for(int i=0; i<N; i++){
      NumericMatrix mat_i = list_of_data[i];
      int Ti = mat_i.nrow();
      arma::mat Xi(Ti, M);
      for(int t=0; t<Ti; t++){
        for(int mm=0; mm<M; mm++){
          Xi(t, mm) = mat_i(t, mm);
        }
      }
      for(int k2=0; k2<K; k2++){
        // diagonal of Psi
        arma::vec psi_diag(M);
        for(int mm=0; mm<M; mm++){
          psi_diag[mm] = Psi[k2](mm,mm);
        }
        // row-wise log-density
        arma::vec lvec = dmvnormFA_rowwiseCpp(Xi, mu[k2], Lambda[k2], psi_diag);
        double sumLog = arma::accu(lvec);
        logLik_i(i, k2) = logPi[k2] + sumLog;
      }
    }

    // Hard assignment
    arma::vec rowMax = arma::max(logLik_i, 1);
    IntegerVector zR_new(N);
    for(int i=0; i<N; i++){
      double bestVal = -1e15;
      int bestk = 1;
      for(int k2=0; k2<K; k2++){
        double val = logLik_i(i,k2);
        if(val > bestVal){
          bestVal = val;
          bestk   = k2 + 1;
        }
      }
      zR_new[i] = bestk;
    }

    // update pi_k
    for(int k2=0; k2<K; k2++){
      int countk = 0;
      for(int i=0; i<N; i++){
        if(zR_new[i] == (k2+1)){
          countk++;
        }
      }
      pi_k[k2] = double(countk) / double(N);
    }

    // log-likelihood
    double new_loglik = arma::accu(rowMax);
    double diff = std::fabs(new_loglik - old_loglik);

    // Check convergence
    // 1) If assignments haven't changed, stop
    bool sameAll = true;
    for(int i=0; i<N; i++){
      if(zR[i] != zR_new[i]){
        sameAll = false;
        break;
      }
    }

    if(iterEM >= 5){
      if(sameAll){
        Rcpp::Rcout << "Converged by assignment at iter=" << iterEM << "\n";
        zR = zR_new;
        break;
      }

      // 2) If log-likelihood diff < tol, stop
      if(diff < tol){
        Rcpp::Rcout << "Converged by loglik diff at iter=" << iterEM
                    << " diff=" << diff << "\n";
        zR = zR_new;
        break;
      }
    }
    zR = zR_new;
    old_loglik = new_loglik;
  }

  // Prepare output
  Rcpp::IntegerVector z_out(N);
  for(int i=0; i<N; i++){
    z_out[i] = zR[i];
  }

  Rcpp::List Lmbd(K), Mu(K), PsiList(K);
  for(int k=0; k<K; k++){
    Lmbd[k]    = wrap(Lambda[k]);
    Mu[k]      = wrap(mu[k]);
    PsiList[k] = wrap(Psi[k]);
  }

  NumericVector pi_out(K);
  for(int k2=0; k2<K; k2++){
    pi_out[k2] = pi_k[k2];
  }

  Rcpp::List ret;
  ret["z"]      = z_out;
  ret["Lambda"] = Lmbd;
  ret["mu"]     = Mu;
  ret["Psi"]    = PsiList;
  ret["pi"]     = pi_out;

  return ret;
}
