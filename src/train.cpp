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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]


using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */


double logsumexp(double &x, double &y, bool flg){
  if (flg) return y; // init mode
  if (x == y) return x + 0.69314718055; // log(2)
  double vmin = std::min(x, y);
  double vmax = std::max(x, y);
  if (vmax > vmin + 50){
    return vmax;
  } else {
    return vmax + std::log(std::exp(vmin - vmax) + 1.0);
  }
}

double logsumexp_Eigen(VectorXd &vec){
  double sum = 0.0;
  for(int i = 0; i < vec.size(); i++){
    sum = logsumexp(sum, vec[i], (i == 0));
  }
  return sum;
}

double gammapdfln(double x, double a, double b){
  return a * log(b) - lgamma(a) + (a - 1.0) * log(x) - b * x;
}

NumericVector alpha_reformat(VectorXd& alpha, int& num_topics){
	NumericVector alpha_rvec(num_topics);

	for(int i=0; i<num_topics; ++i){
		alpha_rvec[i] = alpha(i);
	}

	return alpha_rvec;
}

// int rcat(Eigen::VectorXd &prob){ // was called 'multi1'
//   double u = R::runif(0, 1);
//   double temp = 0.0;
//   int index = 0;
//   for (int ii = 0; ii < prob.size(); ii++){
//     temp += prob(ii);
//     if (u < temp){
//       index = ii;
//       break;
//     }
//   }
//   return index;
// }

// int rcat_without_normalize(Eigen::VectorXd &prob, double &total){ // was called 'multi1'
//   double u = R::runif(0, 1) * total;
//   double temp = 0.0;
//   int index = 0;
//   for (int ii = 0; ii < prob.size(); ii++){
//     temp += prob(ii);
//     if (u < temp){
//       index = ii;
//       break;
//     }
//   }
//   return index;
// }

double loglikelihood_normal(MatrixXi& n_x0_kv,
                      MatrixXi& n_x1_kv,
                      VectorXi& n_x0_k, VectorXi& n_x1_k,
                      MatrixXd& n_dk, VectorXd& alpha,
                      double beta, double beta_s,
                      double gamma_1, double gamma_2,
                      int num_topics, int k_seeded, int num_vocab, int num_doc,
                      std::vector< std::unordered_set<int> > &keywords) {
  double loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += lgamma(beta + (double)n_x0_kv(k, v) ) - lgamma(beta);
      loglik += lgamma(beta_s + (double)n_x1_kv(k, v) ) - lgamma(beta_s);
    }
    // word normalization
    loglik += lgamma( beta * (double)num_vocab ) - lgamma(beta * (double)num_vocab + (double)n_x0_k(k) );
    loglik += lgamma( beta_s * (double)num_vocab ) - lgamma(beta_s * (double)num_vocab + (double)n_x1_k(k) );
    // x
    loglik += lgamma( (double)n_x0_k(k) + gamma_2 ) - lgamma((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2)
      + lgamma( (double)n_x1_k(k) + gamma_1 ) ;

		// Rcout << (double)n_x0_k(k) << " / " << (double)n_x1_k(k) << std::endl; // debug

    // x normalization
    loglik += lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) - lgamma(gamma_2);
  }
  // z
  for (int d = 0; d < num_doc; d++){
    loglik += lgamma( alpha.sum() ) - lgamma( n_dk.row(d).sum() + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += lgamma( n_dk(d,k) + alpha(k) ) - lgamma( alpha(k) );
    }
  }
  return loglik;
}

// takes the sufficient statistics and parameters of the model and
// returns a new value for z
int sample_z(MatrixXi& n_x0_kv,
             MatrixXi& n_x1_kv,
             VectorXi& n_x0_k, VectorXi& n_x1_k,
             MatrixXd& n_dk, VectorXd& alpha,
             std::vector<int>& seed_num,
             int x, int z, int w, int doc_id, int num_vocab,
             int num_topics, int k_seeded,
             double gamma_1, double gamma_2, double beta, double beta_s,
             std::vector< std::unordered_set<int> > &keywords){
  // remove data
  if (x == 0){
    n_x0_kv(z, w) -= 1;
    n_x0_k(z) -= 1;
  } else if (x==1) {
    n_x1_kv(z, w) -= 1;
    n_x1_k(z) -= 1;
  } else {
		Rcerr << "Error at sample_z, remove" << std::endl;
	}

  n_dk(doc_id, z) -= 1;

  VectorXd z_prob_vec = VectorXd::Zero(num_topics);
  int new_z = -1; // debug
  double numerator, denominator;
  if (x == 0){
    for (int k = 0; k < num_topics; ++k){
      numerator = (beta + (double)n_x0_kv(k, w)) *
        ((double)n_x0_k(k) + gamma_2) *
        ((double)n_dk(doc_id, k) + alpha(k));

      denominator = ((double)num_vocab * beta + (double)n_x0_k(k)) *
        ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);

      z_prob_vec(k) = numerator / denominator;
    }

    double sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum); // take a sample

  } else {
    for (int k = 0; k < num_topics; ++k){
      if (keywords[k].find(w) == keywords[k].end()){
        z_prob_vec(k) = 0.0;
				continue;
      } else{ // w not one of the seeds
        numerator = (beta_s + (double)n_x1_kv(k, w)) *
          ( ((double)n_x1_k(k) + gamma_1) ) *
          ( ((double)n_dk(doc_id, k) + alpha(k)) );
      }
      denominator = ((double)seed_num[k] * beta_s + (double)n_x1_k(k) ) *
        ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);

      z_prob_vec(k) = numerator / denominator;
    }


		double sum = z_prob_vec.sum();
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum); // take a sample

  }

  // add back data counts
  if (x == 0){
    n_x0_kv(new_z, w) += 1;
    n_x0_k(new_z) += 1;
  } else if (x==1) {
    n_x1_kv(new_z, w) += 1;
    n_x1_k(new_z) += 1;
  } else {
		Rcerr << "Error at sample_z, add" << std::endl;
	}
  n_dk(doc_id, new_z) += 1;

  return new_z;
}


// takes the sufficient statistics and parameters of the model and
// returns a new value for x
int sample_x(MatrixXi& n_x0_kv,
             MatrixXi& n_x1_kv,
             VectorXi& n_x0_k, VectorXi& n_x1_k,
             MatrixXd& n_dk, VectorXd& alpha,
             std::vector<int>& seed_num,
             int x, int z, int w, int doc_id, int num_vocab,
             int num_topics, int k_seeded,
             double gamma_1, double gamma_2, double beta, double beta_s,
             std::vector< std::unordered_set<int> > &keywords){
  // remove data
  if (x == 0){
    n_x0_kv(z, w) -= 1;
    n_x0_k(z) -= 1;
  } else {
    n_x1_kv(z, w) -= 1;
    n_x1_k(z) -= 1;
  }
  n_dk(doc_id, z) -= 1; // not necessary to remove

  // newprob_x1()
  double x1_prob;
  int k = z;
  double numerator;
  double denominator;
  if ( keywords[k].find(w) == keywords[k].end() ){
       x1_prob = -1.0;
  } else {
    numerator = (beta_s + (double)n_x1_kv(k, w)) *
      ( ((double)n_x1_k(k) + gamma_1) );
    denominator = ((double)seed_num[k] * beta_s + (double)n_x1_k(k) ) *
      ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
    x1_prob = numerator / denominator;
  }

  int new_x;
  if(x1_prob == -1.0){
    // if probability of x_di = 1 case is 0, it should be x=0 (regular topic)
    new_x = 0;
  } else {
    // newprob_x0()
    int k = z;
    numerator = (beta + (double)n_x0_kv(k, w)) *
      ((double)n_x0_k(k) + gamma_2);

    denominator = ((double)num_vocab * beta + (double)n_x0_k(k) ) *
      ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
    double x0_prob = numerator / denominator;

    // Normalize
    double sum = x0_prob + x1_prob;

    x1_prob = x1_prob / sum;
    new_x = R::runif(0,1) <= x1_prob;  //new_x = Bern(x0_prob, x1_prob);
  }
  // add back data counts
  if (new_x == 0){
    n_x0_kv(z, w) += 1;
    n_x0_k(z) += 1;
  } else {
    n_x1_kv(z, w) += 1;
    n_x1_k(z) += 1;
  }
  n_dk(doc_id, z) += 1;

  return new_x;
}


// Wrapper around R's uniform random number generator for std shuffles
// inline int rand_wrapper(const int n) { return floor(unif_rand() * n); }

// Turns K into random permutation of 0 to K-1 using R's RNG
std::vector<int> shuffled_indexes(int m) {
  std::vector<int> v(m);
  std::iota(v.begin(), v.end(), 0);
  std::random_shuffle(v.begin(), v.end(), sampler::rand_wrapper);
  return v;
}

// a version of randgen::uniform(lower, upper) that uses R's RNG
// double slice_uniform(double lower, double upper){
//   return lower + (upper - lower) * unif_rand();
// }

// This used to differentiate between input_alpha and alpha
// for reasons I could not understand.  TODO: check this is still correct?
double alpha_loglik(VectorXd& alpha, MatrixXd& n_dk,
                    int num_topics, int k_seeded, int num_doc,
                    double eta_1, double eta_2, double eta_1_regular, double eta_2_regular){
  double loglik = 0.0;
  double fixed_part = 0.0;
	// double eta_1 = 1;
	// double eta_2 = 1;
	// double eta_1_regular = 2;
	// double eta_2_regular = 1;

  MatrixXd ndk_a = n_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting


  fixed_part += lgamma(alpha.sum()); // first term numerator
  for(int k = 0; k < num_topics; k++){
    fixed_part -= lgamma(alpha(k)); // first term denominator
    // Add prior
		if(k < k_seeded){
			loglik += gammapdfln(alpha(k), eta_1, eta_2);
		}else{
			loglik += gammapdfln(alpha(k), eta_1_regular, eta_2_regular);
		}

  }
  for(int d = 0; d < num_doc; d++){
    loglik += fixed_part;
    // second term numerator
    for(int k = 0; k < num_topics; k++){
      loglik += lgamma(ndk_a(d,k));
    }
    // second term denominator
    loglik -= lgamma(ndk_a.row(d).sum());

  }
  return loglik;
}

// updates alpha in place, currently just hands back alpha (by reference)
void slice_sample_alpha(VectorXd& alpha, MatrixXd& n_dk,
                             int num_topics, int k_seeded, int num_doc,
                             double eta_1, double eta_2, double eta_1_regular,
                             double eta_2_regular,
                             double min_v = 1e-9, double max_v = 100.0,
                             int max_shrink_time = 1000){

  double start, end, previous_p, new_p, newlikelihood, slice_;
  VectorXd keep_current_param = alpha;
  std::vector<int> topic_ids = shuffled_indexes(num_topics);
	double store_loglik = alpha_loglik(alpha, n_dk, num_topics, k_seeded, num_doc, eta_1, eta_2, eta_1_regular, eta_2_regular);
	double newalphallk = 0.0;

  for(int i = 0; i < num_topics; i++){
    int k = topic_ids[i];
    start = min_v / (1.0 + min_v); // shrinkp
    end = 1.0;
    // end = shrinkp(max_v);
		previous_p = alpha(k) / (1.0 + alpha(k)); // shrinkp
    slice_ = store_loglik - 2.0 * log(1.0 - previous_p) 
						+ log(unif_rand()); // <-- using R random uniform

    for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++){
      new_p = sampler::slice_uniform(start, end); // <-- using R function above
      alpha(k) = new_p / (1.0 - new_p); // expandp

			newalphallk = alpha_loglik(alpha, n_dk, num_topics, k_seeded, num_doc, eta_1, eta_2, eta_1_regular, eta_2_regular);
      newlikelihood = newalphallk - 2.0 * log(1.0 - new_p);

      if (slice_ < newlikelihood){
				store_loglik = newalphallk;
        break;
      } else if (previous_p < new_p){
        end = new_p;
      } else if (new_p < previous_p){
        start = new_p;
      } else {
				Rcerr << "Something goes wrong in slice_sample_alpha()" << std::endl;
        alpha(k) = keep_current_param(k);
				break;
      }
    }
  }
}

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


// Sample Lambda
double likelihood_lambda(MatrixXd &Lambda, MatrixXd &C, MatrixXd &n_dk,
		double &mu, double &sigma,
		int &num_doc, int &num_topics, int& num_cov
		)
{
	double loglik = 0.0;
	MatrixXd Alpha = (C * Lambda.transpose()).array().exp();
	VectorXd alpha = VectorXd::Zero(num_topics);

	for(int d=0; d<num_doc; d++){
		alpha = Alpha.row(d).transpose(); // Doc alpha, column vector
		// alpha = ((C.row(d) * Lambda)).array().exp(); // Doc alpha, column vector

		loglik += lgamma(alpha.sum()); 
				// the first term numerator in the first square bracket
		loglik -= lgamma( n_dk.row(d).sum() + alpha.sum() ); 
				// the second term denoinator in the first square bracket

		for(int k=0; k<num_topics; k++){
			loglik -= lgamma(alpha(k));
				// the first term denominator in the first square bracket
			loglik += lgamma( n_dk(d, k) + alpha(k) );
				// the second term numerator in the firist square bracket
		}
	}

	// Prior
	double prior_fixedterm = -0.5 * log(2.0 * PI_V * std::pow(sigma, 2.0) );
	for(int k=0; k<num_topics; k++){
		for(int t=0; t<num_cov; t++){
			loglik += prior_fixedterm;
			loglik -= ( std::pow( (Lambda(k,t) - mu) , 2.0) / (2.0 * std::pow(sigma, 2.0)) );
		}
	}

	return loglik;

}

double expand(double& p, double& A){
	double res = -(1.0/A) * log((1.0/p) - 1.0);
	return res;
}

double shrink(double& x, double& A){
	double res = 1.0 / (1.0 + exp(-A*x));
	return res;
}

void sample_lambda_slice(MatrixXd& Lambda, MatrixXd& C,
									 MatrixXd& n_dk,
                   int& num_topics, int& num_cov,
									 int& k_seeded, int& num_doc,
                   double& mu, double& sigma,
                   int max_shrink_time = 1000)
{
	// Slice sampling for Lambda

	double start, end, previous_p, new_p, newlikelihood, slice_, current_lambda;
  std::vector<int> topic_ids = shuffled_indexes(num_topics);
	std::vector<int> cov_ids = shuffled_indexes(num_cov);
	static double A = 1.2; // important, 0.5-1.5, if the inference goes wrong, change it

	double store_loglik = likelihood_lambda(Lambda, C, n_dk,
			mu, sigma, num_doc, num_topics, num_cov);
	double newlambdallk = 0.0;

	for(int kk=0; kk<num_topics; kk++){
		int k = topic_ids[kk];

		for(int tt=0; tt<num_cov; tt++){
			int t = cov_ids[tt];

			start = 0.0; // shrink
			end = 1.0; // shrink

			current_lambda = Lambda(k,t);
			previous_p = shrink(current_lambda, A);
			slice_ = store_loglik - std::log(A * previous_p * (1.0 - previous_p)) 
							+ log(unif_rand()); // <-- using R random uniform


			for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++){
				new_p = sampler::slice_uniform(start, end); // <-- using R function above
				Lambda(k,t) = expand(new_p, A); // expand

				newlambdallk = likelihood_lambda(Lambda, C, n_dk,
						mu, sigma, num_doc, num_topics, num_cov);

				newlikelihood = newlambdallk - std::log(A * new_p * (1.0 - new_p));

				if (slice_ < newlikelihood){
					store_loglik = newlambdallk;
					break;
				} else if (previous_p < new_p){
					end = new_p;
				} else if (new_p < previous_p){
					start = new_p;
				} else {
					Rcerr << "Something goes wrong in sample_lambda_slice()" << std::endl;
					Lambda(k,t) = current_lambda;
					break;
				}

			} // for loop for shrink time

		} // for loop for num_cov
	} // for loop for num_topics

}


void proposal_lambda(MatrixXd& Lambda,
		int& k, int& num_cov, double& mh_sigma){

	for(int i=0; i<num_cov; i++){
		Lambda.coeffRef(k, i) += R::rnorm(0.0, mh_sigma);
	}

}


void sample_lambda_mh(MatrixXd& Lambda, MatrixXd& C,
									 MatrixXd& n_dk,
                   int& num_topics, int& num_cov,
									 int& k_seeded, int& num_doc,
									 std::vector<int>& mh_info,
									 double& mu, double& sigma,
									 double mh_sigma=0.05)
{
	
	std::vector<int> topic_ids = shuffled_indexes(num_topics);
	VectorXd Lambda_current;
	double llk_current;
	double llk_proposal;
	double diffllk;
	double r, u;

	for(int kk=0; kk<num_topics; kk++){
		int k = topic_ids[kk];
		mh_info[1] += 1; // how many times we run mh

		Lambda_current = Lambda.row(k).transpose();
		
		// Current llk
		llk_current = likelihood_lambda(Lambda, C, n_dk, mu, sigma,
				num_doc, num_topics, num_cov);
		
		// Proposal
		proposal_lambda(Lambda, k, num_cov, mh_sigma);
		llk_proposal = likelihood_lambda(Lambda, C, n_dk, mu, sigma,
				num_doc, num_topics, num_cov);
		
		diffllk = llk_proposal - llk_current;
		r = std::min(0.0, diffllk);
		u = log(unif_rand());
		
		if (u < r){
			mh_info[0] += 1; // number of times accepted	
		}else{
			// Put back original values
			for(int i=0; i<num_cov; i++){
				Lambda.coeffRef(k, i) = Lambda_current(i);
			}
		}
	
	}

}

void sample_lambda_mh_single(MatrixXd& Lambda, MatrixXd& C,
									 MatrixXd& n_dk,
                   int& num_topics, int& num_cov,
									 int& k_seeded, int& num_doc,
									 std::vector<int>& mh_info,
									 double& mu, double& sigma,
									 double mh_sigma=0.05)
{
	
	std::vector<int> topic_ids = shuffled_indexes(num_topics);
	std::vector<int> cov_ids = shuffled_indexes(num_cov);
	double Lambda_current;
	double llk_current;
	double llk_proposal;
	double diffllk;
	double r, u;

	for(int kk=0; kk<num_topics; kk++){
		int k = topic_ids[kk];
		
		for(int tt=0; tt<num_cov; tt++){
			int t = cov_ids[tt];
		
			mh_info[1] += 1; // how many times we run mh

			Lambda_current = Lambda(k,t);
			
			// Current llk
			llk_current = likelihood_lambda(Lambda, C, n_dk, mu, sigma,
					num_doc, num_topics, num_cov);
			
			// Proposal
			Lambda(k, t) += R::rnorm(0.0, mh_sigma);
			llk_proposal = likelihood_lambda(Lambda, C, n_dk, mu, sigma,
					num_doc, num_topics, num_cov);
			
			diffllk = llk_proposal - llk_current;
			r = std::min(0.0, diffllk);
			u = log(unif_rand());
			
			if (u < r){
				mh_info[0] += 1; // number of times accepted	
			}else{
				// Put back original values
				Lambda(k, t) = Lambda_current;
				
			}

		}

	}

}


void sample_lambda(MatrixXd& Lambda, MatrixXd& C,
									 MatrixXd& n_dk,
                   int& num_topics, int& num_cov,
									 int& k_seeded, int& num_doc,
									 std::vector<int>& mh_info,
                   double mu=0.0, double sigma=50.0,  // for prior
                   double min_v = 1e-9, double max_v = 100.0,
                   int max_shrink_time = 1000)
{

	double u = unif_rand(); // select sampling methods randomly

	// if(u < 0.3){
	// 	sample_lambda_mh(Lambda, C, n_dk,
 //                   num_topics, num_cov,
	// 								 k_seeded, num_doc, mh_info, mu, sigma);
	// }else if (u < 0.6){
	// 	sample_lambda_mh_single(Lambda, C, n_dk,
 //                   num_topics, num_cov,
	// 								 k_seeded, num_doc, mh_info, mu, sigma);
	// }else{
	// 	sample_lambda_slice(Lambda, C, n_dk, num_topics, num_cov,
	// 				k_seeded, num_doc, mu, sigma);	
	// }

	if(u < 0.4){
		sample_lambda_mh_single(Lambda, C, n_dk,
                   num_topics, num_cov,
									 k_seeded, num_doc, mh_info, mu, sigma);
	}else{
		sample_lambda_slice(Lambda, C, n_dk, num_topics, num_cov,
					k_seeded, num_doc, mu, sigma);
	}

	// sample_lambda_slice(Lambda, C, n_dk, num_topics, num_cov,
				// k_seeded, num_doc, mu, sigma);
	
}


double loglikelihood_cov(MatrixXi& n_x0_kv,
                      MatrixXi& n_x1_kv,
                      VectorXi& n_x0_k, VectorXi& n_x1_k,
                      MatrixXd& n_dk,
											MatrixXd& Lambda, MatrixXd& C,
                      double beta, double beta_s,
                      double gamma_1, double gamma_2,
                      int num_topics, int k_seeded, int num_vocab, int num_doc,
                      std::vector< std::unordered_set<int> > &keywords) {

  double loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += lgamma(beta + (double)n_x0_kv(k, v) ) - lgamma(beta);
      loglik += lgamma(beta_s + (double)n_x1_kv(k, v) ) - lgamma(beta_s);
    }
    // word normalization
    loglik += lgamma( beta * (double)num_vocab ) - lgamma(beta * (double)num_vocab + (double)n_x0_k(k) );
    loglik += lgamma( beta_s * (double)num_vocab ) - lgamma(beta_s * (double)num_vocab + (double)n_x1_k(k) );
    // x
    loglik += lgamma( (double)n_x0_k(k) + gamma_2 ) - lgamma((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2)
      + lgamma( (double)n_x1_k(k) + gamma_1 ) ;

    // x normalization
    loglik += lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) - lgamma(gamma_2);
  }


  // z
	MatrixXd Alpha = (C * Lambda.transpose()).array().exp();
	VectorXd alpha = VectorXd::Zero(num_topics);

  for (int d = 0; d < num_doc; d++){
		alpha = Alpha.row(d).transpose(); // Doc alpha, column vector	

    loglik += lgamma( alpha.sum() ) - lgamma( n_dk.row(d).sum() + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += lgamma( n_dk(d,k) + alpha(k) ) - lgamma( alpha(k) );
    }
  }
  return loglik;
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

	auto start = std::chrono::high_resolution_clock::now(); // track time
	// Simulation Setting
	int simulation_treatment = 0;

	// Data
  List W = model["W"], Z = model["Z"], X = model["X"];
  StringVector files = model["files"], vocab = model["vocab"];
  NumericVector nv_alpha = model["alpha"];
  double gamma_1 = model["gamma_1"], gamma_2 = model["gamma_2"];
  double beta = model["beta"], beta_s = model["beta_s"];
  int k_free = model["extra_k"];
  List seeds = model["seeds"];
  int k_seeded = seeds.size();
	List model_fit = model["model_fit"];

  // document-related constants
  int num_vocab = vocab.size(), num_doc = files.size();

  // Alpha
  int num_topics = k_seeded + k_free;
	MatrixXd Alpha = MatrixXd::Zero(num_doc, num_topics);
	VectorXd alpha = VectorXd::Zero(num_topics); // use in iteration
	// double alpha_initial_value = 50.0 / num_topics;
	// VectorXd alpha_initial = VectorXd::Constant(num_topics, alpha_initial_value);

	// Covariate
	NumericMatrix C_r = model["C"];
	MatrixXd C = Rcpp::as<Eigen::MatrixXd>(C_r);
	int num_cov = C.cols();

	// Covariate for Experiment (Just for simulation, comment out later)
	// MatrixXd C_treated;
	// MatrixXd C_control;
	// if(simulation_treatment){
	// 	C_treated = MatrixXd::Zero(num_doc, num_cov);
	// 	C_control = MatrixXd::Zero(num_doc, num_cov);
	// 	for(int d=0; d<num_doc; d++){
	// 		for(int c=0; c<num_cov; c++){
	// 			C_treated(d, c) = C(d,c);
	// 			C_control(d, c) = C(d,c);
	//
	// 			if(c==0){
	// 				C_treated(d, c) = 1.0;
	// 				C_control(d, c) = 0.0;
	// 			}
	// 		}	
	// 	}
	// }


	// For causal experiment simulation
	// MatrixXd diff_mean_store = MatrixXd::Zero(iter, num_topics);
	// MatrixXd Alpha_treated = MatrixXd::Zero(num_doc, num_topics);
	// MatrixXd Alpha_control = MatrixXd::Zero(num_doc, num_topics);

	// Lambda
	MatrixXd Lambda = MatrixXd::Zero(num_topics, num_cov);
	for(int k=0; k<num_topics; k++){
		// Initialize with R random
		for(int i=0; i<num_cov; i++){
			Lambda.coeffRef(k, i) = R::rnorm(0.0, 0.3);
		}
	}

	// Vector that stores seed words (words in dictionary)
  std::vector< std::unordered_set<int> > keywords(num_topics);
  std::vector<int> seed_num(num_topics);
  for (int ii = 0; ii < k_seeded; ii++){
    IntegerVector wd_ids = seeds[ii];
    seed_num[ii] = wd_ids.size();
    std::unordered_set<int> keywords_set;
    for (int jj = 0; jj < wd_ids.size(); jj++)
			keywords_set.insert(wd_ids(jj));

    keywords[ii] = keywords_set;
  }
	for(int i=k_seeded; i<num_topics; i++){
		std::unordered_set<int> keywords_set{ -1 };
		seed_num[i] = 0;
		keywords[i] = keywords_set;
	}


  // storage for sufficient statistics and their margins
  MatrixXi n_x0_kv = MatrixXi::Zero(num_topics, num_vocab);
  MatrixXi n_x1_kv = MatrixXi::Zero(num_topics, num_vocab);
  MatrixXd n_dk = MatrixXd::Zero(num_doc, num_topics);
  VectorXi n_x0_k = VectorXi::Zero(num_topics);
  VectorXi n_x1_k = VectorXi::Zero(num_topics);


	// Sampling Information
	std::vector<int> mh_info{0,0};

  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    IntegerVector doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
    for(int w_position = 0; w_position < doc_x.size(); w_position++){
      int x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];
      if (x == 0){
        n_x0_kv(z, w) += 1;
        n_x0_k(z) += 1;
      } else {
        n_x1_kv(z, w) += 1;
        n_x1_k(z) += 1;
      }
      n_dk(doc_id, z) += 1.0;
    }
  }
  int total_words = (int)n_dk.sum();

	double prepare_data = std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::high_resolution_clock::now() - start).count();  // track time


  // Randomized update sequence
  for (int it = 0; it < iter; it++){
    std::vector<int> doc_indexes = shuffled_indexes(num_doc); // shuffle

		// Create Alpha for this iteration
		Alpha = (C * Lambda.transpose()).array().exp();
		
    for (int ii = 0; ii < num_doc; ii++){
      int doc_id = doc_indexes[ii];
      IntegerVector doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
			
      std::vector<int> token_indexes = shuffled_indexes(doc_x.size()); //shuffle
			
			// Prepare Alpha for the doc
			alpha = Alpha.row(doc_id).transpose(); // take out alpha
			
			// Iterate each word in the document
      for (int jj = 0; jj < doc_x.size(); jj++){
        int w_position = token_indexes[jj];
        int x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];
			
        doc_z[w_position] = sample_z(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk,
                                     alpha, seed_num, x, z, w, doc_id,
                                     num_vocab, num_topics, k_seeded, gamma_1,
                                     gamma_2, beta, beta_s, keywords);
			
				z = doc_z[w_position]; // use updated z
			
        doc_x[w_position] = sample_x(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk,
                                    alpha, seed_num, x, z, w, doc_id,
                                    num_vocab, num_topics, k_seeded, gamma_1,
                                    gamma_2, beta, beta_s, keywords);
			
      }
			
			X[doc_id] = doc_x; // is doc_x not a pointer/ref to X[doc_id]?
			Z[doc_id] = doc_z;
			
    }
		sample_lambda(Lambda, C, n_dk, num_topics, num_cov,
				k_seeded, num_doc, mh_info);

		
		// Store Lambda
		Rcpp::NumericMatrix Lambda_R = Rcpp::wrap(Lambda);
		List Lambda_iter = model["Lambda_iter"];
		Lambda_iter.push_back(Lambda_R);
		model["Lambda_iter"] = Lambda_iter;

		// if(simulation_treatment){
		// 	// Treatment Effect (Just for simulation, comment out later)
		// 		// Store Treatment Effect
		// 	Alpha_treated = (C_treated * Lambda.transpose()).array().exp();
		// 	VectorXd E_treated = VectorXd::Zero(num_doc);
		//
		// 	Alpha_control = (C_control * Lambda.transpose()).array().exp();
		// 	VectorXd E_control = VectorXd::Zero(num_doc);
		//
		// 	for(int k=0; k<num_topics; k++){
		// 		E_treated = Alpha_treated.col(k).array() / Alpha_treated.rowwise().sum().array();	
		// 		E_control = Alpha_control.col(k).array() / Alpha_control.rowwise().sum().array();
		// 		VectorXd Diff = E_treated.array() - E_control.array();
		// 		double diff_mean = Diff.mean();
		// 		diff_mean_store.coeffRef(it, k) = diff_mean;
		// 	}
		// }


		// Log-likelihood and Perplexity
		int r_index = it + 1;
		if(r_index % output_per == 0 || r_index == 1 || r_index == iter ){
			double loglik = loglikelihood_cov(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk,
																		 Lambda, C,
																		 beta, beta_s, gamma_1, gamma_2,
																		 num_topics, k_seeded, num_vocab, num_doc, keywords);
			double perplexity = exp(-loglik / (double)total_words);
		
			NumericVector model_fit_vec;
			model_fit_vec.push_back(r_index);
			model_fit_vec.push_back(loglik);
			model_fit_vec.push_back(perplexity);
			model_fit.push_back(model_fit_vec);
		
			Rcerr << "[" << r_index << "] log likelihood: " << loglik <<
							 " (perplexity: " << perplexity << ")" << std::endl;
		
		}
		
    checkUserInterrupt();
  }


	model["model_fit"] = model_fit;

	// Add Sampling Info
	Rcpp::IntegerVector sampling_info = Rcpp::wrap(mh_info);
	List sampling_info_list = model["sampling_info"];
	sampling_info_list.push_back(sampling_info);
	model["sampling_info"] = sampling_info_list;

	// Rcpp::NumericMatrix diffmean;
	// if(simulation_treatment){
	// 	// Add Diff mean (comment out later)
	// 	diffmean = Rcpp::wrap(diff_mean_store);
	// 	List sampling_info_list = model["sampling_info"];
	// 	sampling_info_list.push_back(diffmean);
	// 	model["sampling_info"] = sampling_info_list;
	// }

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

