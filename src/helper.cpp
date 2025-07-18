#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;
using namespace std;

//' Calculate the probability for Polya-Gamma Covariate Model
//'
//' Same as utils::calc_PGtheta, but this is for calling from R
//'
//' @param theta_tilda Parameters
//' @param theta Parameters
//' @param num_doc Number of documents
//' @param num_topics Number of topics
//'
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix calc_PGtheta_R(const NumericMatrix &theta_tilda,
                             Eigen::MatrixXd &theta, const int num_doc,
                             const int num_topics) {
  double remaining = 1.0;

  for (int d = 0; d < num_doc; ++d) {
    remaining = 1.0;

    for (int k = 0; k < num_topics; ++k) {
      if (k == 0) {
        theta(d, 0) = theta_tilda(d, 0);
        remaining *= (1 - theta_tilda(d, 0));
      } else if (k == num_topics - 1) {
        theta(d, num_topics - 1) = 1 - theta.row(d).sum();
      } else {
        theta(d, k) = remaining * theta_tilda(d, k);
        remaining *= (1 - theta_tilda(d, k));
      }
    }
  }

  NumericMatrix theta_R = Rcpp::wrap(theta);
  return theta_R;
}
