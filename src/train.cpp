#include <Rcpp.h>

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_set>

// Sampler
#include "sampler.h"

// keyATM models
#include "keyATM_HMM.h"
#include "keyATM_base.h"
#include "keyATM_cov.h"
#include "keyATM_covPG.h"

// Weighted LDA models
#include "LDA_weight.h"
#include "LDA_weightCov.h"
#include "LDA_weightHMM.h"

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;
using namespace std;

//' Run the Collapsed Gibbs sampler for the keyATM Base
//'
//' @param model A initialized model
//' @param resume resume or not
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_base(List model, bool resume = false) {
  keyATMbase keyATMbase_model(model);
  if (resume) {
    keyATMbase_model.resume_fit();
  } else {
    keyATMbase_model.fit();
  }
  model = keyATMbase_model.return_model();
  return model;
}

//' Run the Collapsed Gibbs sampler for the keyATM covariates (Dir-Multi)
//'
//' @param model A initialized model
//' @param resume resume or not
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_cov(List model, bool resume = false) {
  keyATMcov keyATMcov_model(model);
  if (resume) {
    keyATMcov_model.resume_fit();
  } else {
    keyATMcov_model.fit();
  }
  model = keyATMcov_model.return_model();
  return model;
}

//' Run the Collapsed Gibbs sampler for the keyATM covariates (Polya-Gamma)
//'
//' @param model A initialized model
//' @param resume resume or not
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_covPG(List model, bool resume = false) {
  keyATMcovPG keyATMcov_modelPG(model);
  if (resume) {
    Rcout << "Resume is not supported for Polya-Gamma model" << endl;
  } else {
    keyATMcov_modelPG.fit();
  }
  model = keyATMcov_modelPG.return_model();
  return model;
}

//' Run the Collapsed Gibbs sampler for the keyATM Dynamic
//'
//' @param model A initialized model
//' @param resume resume or not
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_HMM(List model, bool resume = false) {
  keyATMhmm hmm_model(model);
  if (resume) {
    hmm_model.resume_fit();
  } else {
    hmm_model.fit();
  }
  model = hmm_model.return_model();
  return model;
}

//' Run the Collapsed Gibbs sampler for weighted LDA
//'
//' @param model A initialized model
//' @param resume resume or not
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_LDA(List model, bool resume = false) {
  LDAweight LDAweight_model(model);
  if (resume) {
    LDAweight_model.resume_fit();
  } else {
    LDAweight_model.fit();
  }
  model = LDAweight_model.return_model();
  return model;
}

//' Run the Collapsed Gibbs sampler for weighted LDA with covariates
//'
//' @param model A initialized model
//' @param resume resume or not
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_LDAcov(List model, bool resume = false) {
  LDAcov ldacov_model(model);
  if (resume) {
    ldacov_model.resume_fit();
  } else {
    ldacov_model.fit();
  }
  model = ldacov_model.return_model();
  return model;
}

//' Run the Collapsed Gibbs sampler for the weighted LDA with HMM model
//'
//' @param model A initialized model
//' @param resume resume or not
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_LDAHMM(List model, bool resume = false) {
  LDAhmm ldahmm_model(model);
  if (resume) {
    ldahmm_model.resume_fit();
  } else {
    ldahmm_model.fit();
  }
  model = ldahmm_model.return_model();
  return model;
}
