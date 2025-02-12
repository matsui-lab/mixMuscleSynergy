// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

/************************************************************
 * dmvnormPCA_rowwiseCpp
 ************************************************************/
// [[Rcpp::export]]
arma::vec dmvnormPCA_rowwiseCpp(const arma::mat &X,
                                const arma::vec &mu,
                                const arma::mat &W,
                                double sigma2)
{
  int T = X.n_rows;
  int M = X.n_cols;

  // Sigma = W W^T + sigma2 I
  mat Sigma = W * W.t();
  Sigma.diag() += sigma2;
  // 数値安定
  Sigma.diag() += 1e-8;

  mat L = chol(Sigma, "lower");

  double logdetSigma = 0.0;
  for(int m=0; m<M; m++){
    logdetSigma += 2.0 * std::log(L(m,m));
  }
  double Mlog2pi = M * std::log(2.0 * M_PI);

  // invert L
  mat Linv = inv(trimatl(L));

  arma::vec out(T, fill::zeros);
  for(int i=0; i<T; i++){
    vec diff = X.row(i).t() - mu;
    vec z = Linv * diff;
    double quadform = dot(z,z);
    double val = -0.5 * (Mlog2pi + logdetSigma + quadform);
    out[i] = val;
  }
  return out;
}

/************************************************************
 * pcaEMupdateCpp
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List pcaEMupdateCpp(const arma::mat &X, int r, int nIter)
{
  int n = X.n_rows;
  int M = X.n_cols;

  // mean
  rowvec mu_r = mean(X, 0);
  mat Xc = X;
  for(int i=0; i<n; i++){
    Xc.row(i) -= mu_r;
  }

  // init W by SVD
  mat U, V;
  vec S;
  svd(U, S, V, Xc);
  mat Vr = V.cols(0, r-1);
  vec sr = S.subvec(0, r-1);
  mat W = Vr * diagmat(sr / std::sqrt((double)n));

  // init sigma^2
  double sigma2 = 1.0;
  if(r < M){
    double extraVar = 0.0;
    for(int j=r; j<M; j++){
      extraVar += std::pow(S[j], 2.0) / double(n);
    }
    extraVar /= double(M - r);
    sigma2 = extraVar;
  }

  // EM
  for(int it=0; it<nIter; it++){
    mat Sigma = W * W.t();
    Sigma.diag() += sigma2;
    Sigma.diag() += 1e-8;

    mat L = chol(Sigma, "lower");
    mat invL = inv(trimatl(L));
    mat invSigma = invL.t() * invL;  // (M x M)

    mat WtSi = W.t() * invSigma;     // (r x M)
    mat Mmat = eye(r, r) + WtSi * W; // (r x r)
    mat M_chol = chol(Mmat, "lower");
    mat M_inv = inv(trimatu(M_chol));
    M_inv = M_inv * M_inv.t();

    mat EZ = M_inv * WtSi * Xc.t();  // (r x n)
    mat EZZt = (EZ * EZ.t()) / double(n);
    EZZt += M_inv;

    mat XcEZt = (Xc.t() * EZ.t()) / double(n);
    mat EZZt_inv = inv(EZZt + 1e-12 * eye(r,r));
    mat W_new = XcEZt * EZZt_inv;

    // sigma2 update
    mat Sxx = (Xc.t() * Xc) / double(n);
    mat tmp = W_new * XcEZt.t();
    mat residual = Sxx - tmp;
    double traceVal = accu(residual.diag());
    double new_sigma2 = traceVal / double(M);
    if(new_sigma2 < 1e-12) new_sigma2 = 1e-12;

    W = W_new;
    sigma2 = new_sigma2;
  }

  Rcpp::List out;
  out["mu"]     = mu_r.t();
  out["W"]      = W;
  out["sigma2"] = sigma2;
  return out;
}

/************************************************************
 * pcaClosedFormCpp
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List pcaClosedFormCpp(const arma::mat &X, int r)
{
  int n = X.n_rows;
  int M = X.n_cols;

  // mean
  rowvec mu_r = mean(X, 0);
  mat Xc = X;
  for(int i=0; i<n; i++){
    Xc.row(i) -= mu_r;
  }

  // SVD
  mat U, V;
  vec S;
  svd(U, S, V, Xc);

  mat Vr = V.cols(0, r-1);
  vec Sr = S.subvec(0, r-1);

  // leftover var => sigma2
  double sigma2 = 0.0;
  if(r < M){
    double sumLeft = 0.0;
    for(int j=r; j<M; j++){
      double lam_j = (S[j]*S[j]) / double(n);
      sumLeft += lam_j;
    }
    sigma2 = sumLeft / double(M - r);
  } else {
    sigma2 = 1e-12;
  }

  // W
  mat W = arma::zeros(M, r);
  for(int j=0; j<r; j++){
    double lam_j = (Sr[j]*Sr[j]) / double(n);
    double diff_j = lam_j - sigma2;
    if(diff_j < 1e-12) diff_j = 1e-12;
    double scale_j = std::sqrt(diff_j);
    W.col(j) = Vr.col(j) * scale_j;
  }

  Rcpp::List out;
  out["mu"]     = mu_r.t();
  out["W"]      = W;
  out["sigma2"] = sigma2;
  return out;
}

/************************************************************
 * mpcaTimeseriesCpp:
 *   ハード割当Mixture PCA
 *   初期クラスタ z_init があれば使用、なければ (i%K)+1
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List mpcaTimeseriesCpp(Rcpp::List list_of_data,
                             int K,
                             int r,
                             int max_iter,
                             int nIterPCA,
                             double tol,
                             std::string method = "EM",
                             Rcpp::Nullable<Rcpp::IntegerVector> z_init = R_NilValue)
{
  int N = list_of_data.size();
  if(N < 1) stop("No data in list_of_data");

  // get M from the first item
  NumericMatrix firstMat = list_of_data[0];
  int M = firstMat.ncol();

  // initial cluster assignment
  IntegerVector zR(N);
  if(z_init.isNotNull()){
    Rcpp::IntegerVector zz(z_init.get());
    if(int(zz.size()) != N){
      stop("z_init size != N");
    }
    for(int i=0; i<N; i++){
      if(zz[i] < 1 || zz[i] > K){
        stop("z_init has out-of-range cluster label");
      }
      zR[i] = zz[i];
    }
  } else {
    for(int i=0; i<N; i++){
      zR[i] = (i % K) + 1;
    }
  }

  // parameters
  std::vector<arma::mat> W(K);
  std::vector<arma::vec> mu(K);
  std::vector<double> sigma2(K);

  // initialize
  std::srand(123);
  for(int k=0; k<K; k++){
    W[k] = 0.01 * arma::randn<arma::mat>(M, r);
    mu[k]= arma::vec(M, fill::zeros);
    sigma2[k] = 0.1;
  }

  arma::vec pi_k(K, fill::value(1.0 / double(K)));
  double old_loglik = -1e15;

  for(int iterEM = 1; iterEM <= max_iter; iterEM++){

    // ============ M-step ============
    for(int k2=0; k2<K; k2++){
      std::vector<int> idx;
      idx.reserve(N);
      for(int i=0; i<N; i++){
        if(zR[i] == (k2+1)) {
          idx.push_back(i);
        }
      }
      if(idx.empty()) {
        // skip
        continue;
      }

      // concat data
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

      // single-cluster PPCA
      Rcpp::List pcaRes;
      if(method == "EM"){
        pcaRes = pcaEMupdateCpp(X_k, r, nIterPCA);
      } else if(method=="closed_form"){
        pcaRes = pcaClosedFormCpp(X_k, r);
      } else {
        Rcpp::stop("method must be 'EM' or 'closed_form'.");
      }

      mu[k2]     = as<arma::vec>(pcaRes["mu"]);
      W[k2]      = as<arma::mat>(pcaRes["W"]);
      sigma2[k2] = as<double>(pcaRes["sigma2"]);
    }

    // ============ E-step ============
    arma::mat logLik_i(N, K, fill::zeros);
    arma::vec logPi(K, fill::zeros);
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
        arma::vec lvec = dmvnormPCA_rowwiseCpp(Xi, mu[k2], W[k2], sigma2[k2]);
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
        double val = logLik_i(i, k2);
        if(val > bestVal){
          bestVal = val;
          bestk   = k2+1;
        }
      }
      zR_new[i] = bestk;
    }

    // pi_k update
    for(int k2=0; k2<K; k2++){
      int countk=0;
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

    // check for assignment convergence
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

  // final output
  Rcpp::IntegerVector z_out(N);
  for(int i=0; i<N; i++){
    z_out[i] = zR[i];
  }

  Rcpp::List W_out(K), Mu_out(K);
  NumericVector sigma2_out(K), pi_out(K);
  for(int k=0; k<K; k++){
    W_out[k]       = wrap(W[k]);
    Mu_out[k]      = wrap(mu[k]);
    sigma2_out[k]  = sigma2[k];
    pi_out[k]      = pi_k[k];
  }

  Rcpp::List ret;
  ret["z"]      = z_out;
  ret["W"]      = W_out;
  ret["mu"]     = Mu_out;
  ret["sigma2"] = sigma2_out;
  ret["pi"]     = pi_out;

  return ret;
}
