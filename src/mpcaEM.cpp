// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

/************************************************************
 * dmvnormPCA_rowwiseCpp:
 *   Sigma = W*W^T + sigma2 * I  を用いた多次元正規の対数密度
 *   (サンプル行ごと) を計算してベクトル返却（シングルスレッド）
 ************************************************************/
// [[Rcpp::export]]
arma::vec dmvnormPCA_rowwiseCpp(const arma::mat &X,
                                const arma::vec &mu,
                                const arma::mat &W,
                                double sigma2)
{
  int T = X.n_rows;  // サンプル数
  int M = X.n_cols;  // 次元
  
  // 共分散行列 Sigma = W W^T + sigma2 * I
  mat Sigma = W * W.t();
  Sigma.diag() += sigma2;
  // 数値安定化
  Sigma.diag() += 1e-8;
  
  // chol
  mat L = chol(Sigma, "lower");
  
  // logdet(Sigma)
  double logdetSigma = 0.0;
  for(int m=0; m<M; m++){
    logdetSigma += 2.0 * std::log(L(m,m));
  }
  double log2pi = std::log(2.0 * M_PI);
  double Mlog2pi = M * log2pi;
  
  // L^-1
  mat Linv = inv(trimatl(L));
  
  // 出力
  arma::vec out(T, fill::zeros);
  
  for(int i=0; i<T; i++){
    vec diff = X.row(i).t() - mu;  // (M×1)
    vec z = Linv * diff;          // (M×1)
    double quadform = dot(z,z);
    double val = -0.5 * (Mlog2pi + logdetSigma + quadform);
    out[i] = val;
  }
  
  return out;
}

/************************************************************
 * pcaEMupdateCpp:
 *   単一クラスタ PPCA を EM で推定する関数
 *   - 引数:
 *       X: (n×M) 行列データ
 *       r: 潜在次元
 *       nIter: EM反復回数
 *   - 戻り値: List (mu, W, sigma2)
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List pcaEMupdateCpp(const arma::mat &X, int r, int nIter)
{
  int n = X.n_rows;  // サンプル数
  int M = X.n_cols;  // 次元
  
  // 平均ベクトル
  rowvec mu_r = mean(X, 0);  
  // 中心化
  mat Xc = X;
  for(int i=0; i<n; i++){
    Xc.row(i) -= mu_r;
  }
  
  // 初期 W => PCA的に SVD (Xc = U S V^T)
  mat U, V;
  vec S;
  svd(U, S, V, Xc);
  mat Vr = V.cols(0, r-1);            // (M×r)
  vec sr = S.subvec(0, r-1);         // r個
  mat W = Vr * diagmat(sr / std::sqrt((double)n));  // (M×r)
  
  // 初期 sigma^2
  double sigma2 = 1.0;
  if(r < M){
    double extraVar = 0.0;
    for(int j=r; j<M; j++){
      extraVar += std::pow(S[j], 2.0) / double(n);
    }
    extraVar /= double(M - r);
    sigma2 = extraVar;
  }
  
  // EM反復
  for(int it=0; it<nIter; it++){
    // Sigma = W W^T + sigma^2 I
    mat Sigma = W * W.t();
    Sigma.diag() += sigma2;
    // 数値安定
    Sigma.diag() += 1e-8;
    
    mat L = chol(Sigma, "lower");
    mat invL = inv(trimatl(L));
    mat invSigma = invL.t() * invL;  // (M×M)
    
    // Mmat = I_r + W^T * invSigma * W
    mat WtSi = W.t() * invSigma;           // (r×M)
    mat Mmat = eye(r, r) + WtSi * W;       // (r×r)
    mat M_chol = chol(Mmat, "lower");
    mat M_inv = inv(trimatu(M_chol));      
    M_inv = M_inv * M_inv.t();            // Mmat^-1
    
    // E[Z] = Mmat^-1 * W^T * invSigma * Xc^T
    mat EZ = M_inv * WtSi * Xc.t();       // (r×n)
    
    // E[Z Z^T] = (1/n) EZ * EZ^T + M_inv
    mat EZZt = (EZ * EZ.t()) / double(n);
    EZZt += M_inv;
    
    // W_new
    mat XcEZt = (Xc.t() * EZ.t()) / double(n); // (M×r)
    mat EZZt_inv = inv(EZZt + 1e-12 * eye(r,r));
    mat W_new = XcEZt * EZZt_inv;
    
    // sigma2_new
    mat Sxx = (Xc.t() * Xc) / double(n); // (M×M)
    mat tmp = W_new * XcEZt.t();        // (M×M)
    mat residual = Sxx - tmp;
    double traceVal = arma::accu(residual.diag());
    double new_sigma2 = traceVal / double(M);
    if(new_sigma2 < 1e-12) new_sigma2 = 1e-12;
    
    // 更新
    W = W_new;
    sigma2 = new_sigma2;
  }
  
  List out;
  out["mu"]     = mu_r.t();  // (M×1)
  out["W"]      = W;         // (M×r)
  out["sigma2"] = sigma2;    // スカラー
  return out;
}

/************************************************************
 * pcaClosedFormCpp:
 *   単一クラスタ PPCA の "解析解(閉形式)" を
 *   SVD から直接求める関数
 *   - 引数:
 *       X: (n×M)
 *       r: 潜在次元
 *   - 戻り値: List (mu, W, sigma2)
 *
 *   ※ Bishopの論文(Tipping & Bishop 1999) などで知られる
 *     PPCAの解析解を用いた実装例
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List pcaClosedFormCpp(const arma::mat &X, int r)
{
  int n = X.n_rows;
  int M = X.n_cols;
  
  // 平均ベクトル
  rowvec mu_r = mean(X, 0);
  
  // 中心化
  mat Xc = X;
  for(int i=0; i<n; i++){
    Xc.row(i) -= mu_r;
  }
  
  // SVD: Xc = U S V^T
  //  S: singular values (>=0)
  mat U, V;
  vec S;
  svd(U, S, V, Xc);
  
  // 寄与度の大きい r 成分を取り出す
  mat Vr = V.cols(0, r-1);         // (M×r)
  vec Sr = S.subvec(0, r-1);       // 長さ r
  
  // 固有値 λ_j = (S_j^2 / n)
  // 残差の分散 sigma^2 = 平均 (λ_j) for j=r..(M-1)
  double sigma2 = 0.0;
  if(r < M){
    double sumLeft = 0.0;
    for(int j=r; j<M; j++){
      double lam_j = (S[j]*S[j]) / double(n);
      sumLeft += lam_j;
    }
    sigma2 = sumLeft / double(M - r);
  } else {
    // r = M の場合は残差ゼロにする
    sigma2 = 1e-12;  // 数値的に完全0だと不安定になりやすいので微小値
  }
  
  // W = V_r * (Lambda_r - sigma^2 I)^(1/2)
  //    ただし Lambda_r = diag( S_r^2 / n )
  //    要素ごとに sqrt(...) が負にならないよう注意
  mat W = arma::zeros(M, r);
  for(int j=0; j<r; j++){
    double lam_j = (Sr[j]*Sr[j]) / double(n); // S_j^2 / n
    double diff_j = lam_j - sigma2;
    if(diff_j < 1e-12) diff_j = 1e-12; // 数値的補正
    double scale_j = std::sqrt(diff_j);
    // W.col(j) = V_r.col(j) * scale_j
    W.col(j) = Vr.col(j) * scale_j;
  }
  
  // 出力
  List out;
  out["mu"]     = mu_r.t(); // (M×1)
  out["W"]      = W;        // (M×r)
  out["sigma2"] = sigma2;
  return out;
}

/************************************************************
 * mpcaTimeseriesCpp:
 *   複数の時系列(行列)からなる list_of_data に対して
 *   Kクラスタの混合 PPCA を推定
 *
 *   - method = "EM"        => pcaEMupdateCpp を用いる
 *   - method = "closed_form" => pcaClosedFormCpp を用いる
 *
 *   ※ "closed_form" 選択時は、M-stepで各クラスタに属するデータを
 *     一括で縦連結し、SVDから直接パラメータを求めます(1回で終わり)。
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List mpcaTimeseriesCpp(Rcpp::List list_of_data,
                             int K,
                             int r,
                             int max_iter,
                             int nIterPCA,
                             double tol,
                             std::string method = "EM")
{
  int N = list_of_data.size();  // データ数(系列数)
  if(N < 1) stop("No data in list_of_data");
  
  // 最初のデータから列数(次元)を取得
  NumericMatrix firstMat = list_of_data[0];
  int M = firstMat.ncol();
  
  // 初期クラスタ割当
  IntegerVector zR(N);
  for(int i=0; i<N; i++){
    zR[i] = (i % K) + 1; // 1～K
  }
  
  // パラメータ
  std::vector<arma::mat> W(K);
  std::vector<arma::vec> mu(K);
  std::vector<double> sigma2(K);
  
  // ランダム初期
  std::srand(123);
  for(int k=0; k<K; k++){
    W[k] = 0.01 * arma::randn<arma::mat>(M, r); // (M×r)
    mu[k] = arma::vec(M, fill::zeros);
    sigma2[k] = 0.1;
  }
  // クラスタ比
  arma::vec pi_k(K, fill::value(1.0 / double(K)));
  
  double old_loglik = -1e15;
  
  // EMループ
  for(int iterEM = 1; iterEM <= max_iter; iterEM++){
    
    //======================
    // M-step
    //======================
    for(int k=0; k<K; k++){
      // クラスタ k に属するデータ(行列)のインデックス
      std::vector<int> idx;
      idx.reserve(N);
      for(int i=0; i<N; i++){
        if(zR[i] == (k+1)) {
          idx.push_back(i);
        }
      }
      if(idx.empty()) {
        // 該当データなし -> 更新しない
        continue;
      }
      
      // 該当データを縦連結
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
      
      // ---------------------------
      // method で単一クラスタ推定を切り替え
      // ---------------------------
      if(method == "EM") {
        // pcaEMupdateCpp (EM反復)
        Rcpp::List pcaRes = pcaEMupdateCpp(X_k, r, nIterPCA);
        mu[k]     = as<arma::vec>(pcaRes["mu"]);
        W[k]      = as<arma::mat>(pcaRes["W"]);
        sigma2[k] = as<double>(pcaRes["sigma2"]);
        
      } else if(method == "closed_form") {
        // pcaClosedFormCpp (解析解)
        Rcpp::List pcaRes = pcaClosedFormCpp(X_k, r);
        mu[k]     = as<arma::vec>(pcaRes["mu"]);
        W[k]      = as<arma::mat>(pcaRes["W"]);
        sigma2[k] = as<double>(pcaRes["sigma2"]);
        
      } else {
        Rcpp::stop("method must be either 'EM' or 'closed_form'");
      }
    }
    
    //======================
    // E-step
    //======================
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
        // 対数尤度 (dmvnormPCA_rowwiseCpp)
        arma::vec lvec = dmvnormPCA_rowwiseCpp(Xi, mu[k2], W[k2], sigma2[k2]);
        double sumLog = arma::accu(lvec);
        logLik_i(i,k2) = logPi[k2] + sumLog;
      }
    }
    
    // z 更新 (ハード割当)
    arma::vec rowMax = arma::max(logLik_i, 1); // (N)
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
    
    // pi_k 更新
    for(int k2=0; k2<K; k2++){
      int countk = 0;
      for(int i=0; i<N; i++){
        if(zR_new[i] == (k2+1)) {
          countk++;
        }
      }
      pi_k[k2] = double(countk) / double(N);
    }
    
    // 対数尤度
    double new_loglik = arma::accu(rowMax);
    double diff = std::fabs(new_loglik - old_loglik);
    
    // 収束判定1: 割当が変わらなければ終了
    bool sameAll = true;
    for(int i=0; i<N; i++){
      if(zR[i] != zR_new[i]){
        sameAll = false;
        break;
      }
    }
    if(sameAll){
      Rcpp::Rcout << "Converged by assignment at iter=" << iterEM << "\n";
      zR = zR_new;
      break;
    }
    // 収束判定2: 対数尤度の差分
    if(diff < tol){
      Rcpp::Rcout << "Converged by loglik diff at iter=" << iterEM
                  << " diff=" << diff << "\n";
      zR = zR_new;
      break;
    }
    
    // 更新
    zR = zR_new;
    old_loglik = new_loglik;
  }
  
  // 結果をリスト化
  //   z_out, W_out, mu_out, sigma2_out, pi_out
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
