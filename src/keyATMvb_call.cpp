#include <Rcpp.h>

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>

#include "keyATMvb_main.h"
#include "sampler.h"
#include <algorithm>
#include <iostream>
#include <unordered_set>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;
using namespace std;

//' Run the Variational Bayes for the keyATM models
//'
//' @param model A model
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATMvb_call(List model) {
  keyATMvb keyATMvb_model(model);
  keyATMvb_model.fit();
  model = keyATMvb_model.return_model();
  return model;
}
