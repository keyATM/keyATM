#ifndef __utils__INCLUDED__
#define __utils__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

namespace utils
{
  void calc_PGtheta(const NumericMatrix theta_tilda, MatrixXd &theta, const int num_doc, const int num_topics);
}

#endif

