#include "utils.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

namespace utils
{
  void calc_PGtheta(const NumericMatrix theta_tilda, MatrixXd &theta, const int num_doc, const int num_topics) {
    double remaining = 1.0;
    theta.setZero();

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
  }

}


