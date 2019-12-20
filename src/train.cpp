#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "sampler.h"
#include "keyATM_base.h"
#include "keyATM_cov.h"
#include "keyATM_HMM.h"
#include "LDA_weight.h"
#include "LDA_weightCov.h"
#include "LDA_weightHMM.h"


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]


using namespace Eigen;
using namespace Rcpp;
using namespace std;


//' Run the Collapsed Gibbs sampler for the keyATM Base
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_base(List model, int iter = 0, int output_per = 10)
{
  keyATMbase keyATMbase_model(model, iter, output_per);
  keyATMbase_model.fit();
  model = keyATMbase_model.return_model();
  return model;
}


//' Run the Collapsed Gibbs sampler for the keyATM covariates
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_cov(List model, int iter = 0, int output_per = 10)
{
  keyATMcov keyATMcov_model(model, iter, output_per);
  keyATMcov_model.fit();
  model = keyATMcov_model.return_model();
  return model;
}


//' Run the Collapsed Gibbs sampler for the keyATM Dynamic
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_HMM(List model, int iter = 0, int output_per = 10)
{
  keyATMhmm hmm_model(model, iter, output_per);
  hmm_model.fit();
  model = hmm_model.return_model();
  return model;
}


//' Run the Collapsed Gibbs sampler for weighted LDA
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_LDA(List model, int iter = 0, int output_per = 10)
{
  LDAweight LDAweight_model(model, iter, output_per);
  LDAweight_model.fit();
  model = LDAweight_model.return_model();
  return model;
}


//' Run the Collapsed Gibbs sampler for weighted LDA with covariates
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_LDAcov(List model, int iter = 0, int output_per = 10)
{
  LDAcov ldacov_model(model, iter, output_per);
  ldacov_model.fit();
  model = ldacov_model.return_model();
  return model;
}




//' Run the Collapsed Gibbs sampler for the weighted LDA with HMM model
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @keywords internal
// [[Rcpp::export]]
List keyATM_fit_LDAHMM(List model, int iter = 0, int output_per = 10)
{
  LDAhmm ldahmm_model(model, iter, output_per);
  ldahmm_model.fit();
  model = ldahmm_model.return_model();
  return model;
}



