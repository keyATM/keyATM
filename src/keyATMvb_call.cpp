#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "sampler.h"
#include "keyATMvb_main.h"

// [[Rcpp::plugins(cpp11)]]
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
List keyATMvb_call(List model)
{
  keyATMvb keyATMvb_model(model);
  keyATMvb_model.fit();
  model = keyATMvb_model.return_model();
  return model;
}
