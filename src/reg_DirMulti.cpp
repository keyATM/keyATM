#include <Rcpp.h>
#include <RcppEigen.h>


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]


using namespace Eigen;
using namespace Rcpp;
using namespace std;


double mylgamma(const double x)
{
  // Good approximation when x > 1
  //    x > 1: max abs err: 2.272e-03
  //    x > 0.5: 0.012
  //    x > 0.6: 0.008
  // Abramowitz and Stegun p.257
  
  // return(lgamma(x));  // Debug
  
  if(x < 0.6)
    return (lgamma(x));
  else
    return ((x-0.5)*log(x) - x + 0.91893853320467 + 1/(12*x));
};


double mydigamma(const double x)
{
  // return(R::digamma(x));  // Debug

  if (x <= 3)
    return(R::digamma(x));

  return (log(x) - 1 / (2*x)  // good enough
          // + 1 / (120*pow(x, 4)) 
          // - 1 / (252 * pow(x, 6)) + 1 / (240 * pow(x, 8))
          // - 5 / (660 * pow(x, 10)) + 691 / (32760 * pow(x, 12))
          // - 1 / (12 * pow(x, 14))
         );
}


double mytrigamma(const double x)
{
  // return(R::trigamma(x));  // degbug

  if (x <= 3)
    return(R::trigamma(x));

  return(1 / x + 1 / (2 * pow(x, 2)));

}


//' Dirichlet Multinomial Distribution
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector ddirmnCpp(Eigen::MatrixXd Y, Eigen::MatrixXd Lambda, Eigen::MatrixXd X)
{
  // Original:
  // https://github.com/cran/MGLM/blob/0942d9778d0f230208cf8f04a75956c9fb2c2a3b/R/pdfln.R

  MatrixXd Alpha = (X * Lambda).array().exp();

  double llk = 0.0;
  int num_doc = Y.rows();
  int num_topics = Y.cols();
  VectorXd alpha = VectorXd::Zero(num_topics);
  VectorXd doc_len = Y.rowwise().sum();
  double alpha_sum;

  for (int d = 0; d < num_doc; ++d) {
    alpha = Alpha.row(d).transpose();
    alpha_sum = alpha.sum();

    if (alpha_sum > 1e+08) {
      NumericVector R_llk = NumericVector::create(R_NegInf);
      return(Rcpp::wrap(R_llk));
    } 
    
    llk += mylgamma(doc_len(d) + 1);
    llk += mylgamma(alpha_sum);

    llk -= mylgamma(doc_len(d) + alpha_sum);

    for (int k = 0; k < num_topics; ++k) {
      llk += mylgamma(Y(d, k) + alpha(k)); 

      llk -= mylgamma(Y(d, k) + 1);
      llk -= mylgamma(alpha(k));
    }
  }

  NumericVector R_llk = NumericVector::create(llk);

  return(Rcpp::wrap(R_llk));
}


//' Calculate Hessian
//'
//' @keywords internal
// [[Rcpp::export]]
List objfun_helper(Eigen::MatrixXd Lambda,
                   Eigen::MatrixXd X,
                   Eigen::MatrixXd Y,
                   int d, int p, List Res)
{
  //
  // Prepare
  //
  int num_doc = Y.rows();
  int num_topics = d;

  MatrixXd Alpha = (X * Lambda).array().exp();
  VectorXd m = Y.rowwise().sum();
  VectorXd Alpha_rowsum = Alpha.rowwise().sum();

  // tmpvector and tmpmatrix
  VectorXd tmpvector = Alpha_rowsum + m;
  VectorXd tmpvector2 = Alpha_rowsum;

  MatrixXd tmpmatrix = Alpha.array() + Y.array();
  MatrixXd tmpmatrix2 = Alpha.array() + Y.array();

  for (int d = 0; d < num_doc; d++) {
    tmpvector(d) = mydigamma(tmpvector(d)) - mydigamma(Alpha_rowsum(d));
      // Original checks NaN here

    tmpvector2(d) = mytrigamma(tmpvector2(d)) - mytrigamma(m(d) + Alpha_rowsum(d));

    for (int k = 0; k < num_topics; ++k) {
      tmpmatrix(d, k) = mydigamma(tmpmatrix(d, k)) - mydigamma(Alpha(d, k));

      tmpmatrix2(d, k) = - mytrigamma(tmpmatrix2(d, k)) + mytrigamma(Alpha(d, k));
    }
  }
  tmpmatrix2 = Alpha.array() * tmpmatrix.array() - Alpha.array().pow(2) * tmpmatrix2.array();

  //
  // Hessian
  //

  // Append
  MatrixXd Beta1 = MatrixXd::Zero(Alpha.rows() * p, Alpha.cols());
  for (int p_index = 0; p_index < p; ++p_index) {
    // Block operations: https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
    Beta1.block(p_index * Alpha.rows(), 0, Alpha.rows(), Alpha.cols()) = Alpha; 
  }
  Beta1.resize(X.rows(), p * d);

  MatrixXd x1 = MatrixXd::Zero(X.rows(), X.cols() * d);
  for (int d_index = 0; d_index < d; ++d_index) {
    x1.block(0, d_index * X.cols(), X.rows(), X.cols()) = X; 
  }

  MatrixXd Hessian = Beta1.array() * x1.array();
  MatrixXd tmp = Hessian.array().colwise() * tmpvector2.array();
  Hessian = Hessian.transpose() * tmp;

  int start;
  VectorXd tmp2;
  MatrixXd tmp3;
  for (int i = 1; i <= d; ++i) {
    start = (i - 1) * p + 1;

    tmp2 = Alpha.col(i - 1).array() * tmpvector.array() - tmpmatrix2.col(i - 1).array();
    tmp3 = X.array().colwise() * tmp2.array();

    Hessian.block(start - 1, start - 1, p, p) =
      Hessian.block(start - 1, start - 1, p, p) -
      X.transpose() * tmp3;
  }
  NumericMatrix Res_hessian = Rcpp::wrap(-Hessian);


  //
  // Grad
  //
  MatrixXd tmp4 = tmpmatrix.colwise() - tmpvector;
  MatrixXd dalpha = Alpha.array() * tmp4.array();

  // Append
  MatrixXd dalpha2 = MatrixXd::Zero(dalpha.rows() * p, dalpha.cols());
  for (int p_index = 0; p_index < p; ++p_index) {
    dalpha2.block(p_index * dalpha.rows(), 0, dalpha.rows(), dalpha.cols()) = dalpha; 
  }
  dalpha2.resize(X.rows(), p * d);

  VectorXd dl = (dalpha2.array() * x1.array()).colwise().sum();
  NumericVector Res_dl = Rcpp::wrap(-dl);
  NumericVector Res_tmpvec = Rcpp::wrap(tmpvector);
  NumericMatrix Res_tmpmat = Rcpp::wrap(tmpmatrix);

  Res["Hessian"] = Res_hessian;
  Res["dl"] = Res_dl;
  Res["tmpvec"] = Res_tmpvec;
  Res["tmpmat"] = Res_tmpmat;

  return(Res);
}

