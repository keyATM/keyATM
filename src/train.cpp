#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "lda_cov.h"
#include "idealpoint.h"
#include "sampler.h"
#include "keyATM_basic.h"
#include "keyATM_cov.h"
#include "keyATM_HMM.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]


using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */



//' Run the Collapsed Gibbs sampler for the standard model
//'
//' @param model A model, from \code{init} or a previous invocation of \code{train}
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List topicdict_train(List model, int iter = 0, int output_per = 10){

	keyATMbasic keyATMbasic_model(model, iter, output_per);
	model = keyATMbasic_model.return_model();
	return model;

}


//' Run the Collapsed Gibbs sampler for the covariate model
//'
//' @param model A model, from \code{init} or a previous invocation of \code{train}, including a covariate
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List topicdict_train_cov(List model, int iter = 0, int output_per = 10){

	keyATMcov keyATMcov_model(model, iter, output_per);
	model = keyATMcov_model.return_model();
	return model;

}


//' Run the Collapsed Gibbs sampler for the HMM model
//'
//' @param model A model, from \code{init} or a previous invocation of \code{train}, including a covariate
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List topicdict_train_HMM(List model, int iter = 0, int output_per = 10){

	// keyATMhmm keyATMhmm_model(model, iter, output_per);
	// model = keyATMhmm_model.return_model();
	return model;

}


//' Run the Collapsed Gibbs sampler for LDA Dir-Multi (Mimno and McCalum 2008)
//'
//' @param model A model, from \code{init} or a previous invocation of \code{train}, including a covariate
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List lda_cov(List model, int K, int iter=0, int output_iter=10)
{

	LDACOV ldacov(model, K, iter, output_iter);

	return model;
}


//' Run the Collapsed Gibbs sampler for Ideal Point Estimation Model
//'
//' @param model A model, from \code{init} or a previous invocation of \code{train}, including a covariate
//' @param author_info author information
//' @param iter Required number of iterations
//' @param output_per Show log-likelihood and perplexity per this number during the iteration
//'
//' @export
// [[Rcpp::export]]
List topicdict_idealpoint(List model, List author_info, int iter=0, int output_iter=10)
{

	IDEALPOINT idealpoint(model, author_info, iter, output_iter);

	return model;
}

