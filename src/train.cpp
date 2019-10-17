#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "lda_cov.h"
#include "sampler.h"
#include "keyATM_basic.h"
#include "keyATM_cov.h"
#include "keyATM_HMM.h"
#include "LDA_weight.h"
#include "LDA_weightHMM.h"


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]


using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */


//' Run the Collapsed Gibbs sampler for the standard model
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List keyATM_train(List model, int iter = 0, int output_per = 10){

  keyATMbasic keyATMbasic_model(model, iter, output_per);
  keyATMbasic_model.fit();
  model = keyATMbasic_model.return_model();
  return model;

}


//' Run the Collapsed Gibbs sampler for the covariate model
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List keyATM_train_cov(List model, int iter = 0, int output_per = 10){

  keyATMcov keyATMcov_model(model, iter, output_per);
  keyATMcov_model.fit();
  model = keyATMcov_model.return_model();
  return model;

}


//' Run the Collapsed Gibbs sampler for the HMM model
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List keyATM_train_HMM(List model, int iter = 0, int output_per = 10){

  keyATMhmm hmm_model(model, iter, output_per);
  hmm_model.fit();
  model = hmm_model.return_model();
  return model;

}


//' Run the Collapsed Gibbs sampler for LDA Dir-Multi (Mimno and McCalum 2008)
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List lda_cov(List model, int K, int iter = 0, int output_iter = 10)
{

  LDACOV ldacov(model, K, iter, output_iter);
  // ldacov.fit();

  return model;
}


//' Run the Collapsed Gibbs sampler for LDA with weights
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List LDA_weight(List model, int iter = 0, int output_per = 10){

  LDAweight LDAweight_model(model, iter, output_per);
  LDAweight_model.fit();
  model = LDAweight_model.return_model();
  return model;

}


//' Run the Collapsed Gibbs sampler for the HMM model
//'
//' @param model A initialized model
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List keyATM_train_LDAHMM(List model, int iter = 0, int output_per = 10){

  LDAhmm ldahmm_model(model, iter, output_per);
  ldahmm_model.fit();
  model = ldahmm_model.return_model();
  return model;

}



