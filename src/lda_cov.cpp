#include "lda_cov.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

// namespace sampler{
// 	inline int rand_wrapper(const int n) { return floor(unif_rand() * n); }
//
// 	double slice_uniform(double& lower, double& upper){
// 		return lower + (upper - lower) * unif_rand();
// 	}
// }

LDACOV::LDACOV(List model_, const int K, const int iter_, const int output_iter_)
{
	// Get Info from R
	model = model_;
	W = model["W"];
	Z = model["Z"];
	num_topics = K;
	C = Rcpp::as<Eigen::MatrixXd>(model["C"]);
  gamma_1 = model["gamma_1"];
	gamma_2 = model["gamma_2"];
  beta = model["beta"];
	iter = iter_;
	output_iter = output_iter_;
	
	vocab = model["vocab"];
	num_vocab = vocab.size();
	num_cov = C.cols();
	num_doc = C.rows();
	
	// Initialize
	initialize();
	
	// Iteration
	iteration();
	
}

void LDACOV::initialize()
{
	// Initialize Lambda
	Lambda = MatrixXd::Zero(num_topics, num_cov);
	for(int k=0; k<num_topics; k++){
		// Initialize with R random
		for(int i=0; i<num_cov; i++){
			Lambda(k, i) = R::rnorm(0.0, 1.5);
		}
	}
	
	// Organize docs
  n_kv = MatrixXi::Zero(num_topics, num_vocab);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_k = VectorXi::Zero(num_topics);
		
  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    IntegerVector doc_z = Z[doc_id], doc_w = W[doc_id];
    for(int w_position = 0; w_position < doc_z.size(); w_position++){
      int z = doc_z[w_position], w = doc_w[w_position];
			n_kv(z, w) += 1;
			n_k(z) += 1;
      n_dk(doc_id, z) += 1.0;
    }
  }
  total_words = (int)n_dk.sum();
	
	Alpha = MatrixXd::Zero(num_doc, num_topics);
}


void LDACOV::iteration()
{
	// Iteration
	VectorXd alpha_ = VectorXd::Zero(num_topics);
	
	for(int it=0; it<iter; it++){
		std::vector<int> doc_indexes = shuffled_indexes(num_doc);
	
		Alpha = (C * Lambda.transpose()).array().exp();
	
    for (int ii = 0; ii < num_doc; ii++){
      int doc_id_ = doc_indexes[ii];
      IntegerVector doc_z = Z[doc_id_], doc_w = W[doc_id_];
			
      std::vector<int> token_indexes = shuffled_indexes(doc_z.size()); //shuffle
			
			// Prepare Alpha for the doc
			alpha_ = Alpha.row(doc_id_).transpose(); // take out alpha
			
			// Iterate each word in the document
      for (int jj = 0; jj < doc_z.size(); jj++){
        int w_position = token_indexes[jj];
        int z_ = doc_z[w_position], w_ = doc_w[w_position];
			
				int new_z = sample_z(alpha_, z_, w_, doc_id_);
        doc_z[w_position] = new_z;
			
      }
			
			Z[doc_id_] = doc_z;
    }
		lambda_sample();
		lambda_store();
	
		int r_index = it + 1;
		if(r_index % output_iter == 0 || r_index == 1 || r_index == iter ){
			loglik_store(r_index);
		}


		checkUserInterrupt();
	}
}

double LDACOV::loglik_lambda()
{
	double loglik = 0.0;
	double mu = 0.0;
	double sigma = 50.0;
	Alpha = (C * Lambda.transpose()).array().exp();
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

void LDACOV::lambda_sample()
{

	// Slice sampling for Lambda

	double start, end, previous_p, new_p, newlikelihood, slice_, current_lambda;
  std::vector<int> topic_ids = shuffled_indexes(num_topics);
	std::vector<int> cov_ids = shuffled_indexes(num_cov);
	static double A = 0.8; // important, 0.5-1.5
	static int max_shrink_time = 1000;

	double store_loglik = loglik_lambda();
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

				newlambdallk = loglik_lambda();

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

	model["model_fit"] = model_fit;

}

void LDACOV::loglik_store(int& r_index)
{
	double loglik = loglik_calc();
	double perplexity = exp(-loglik / (double)total_words);

	NumericVector model_fit_vec;
	model_fit_vec.push_back(r_index);
	model_fit_vec.push_back(loglik);
	model_fit_vec.push_back(perplexity);
	model_fit.push_back(model_fit_vec);

	Rcerr << "[" << r_index << "] log likelihood: " << loglik <<
					 " (perplexity: " << perplexity << ")" << std::endl;
}

double LDACOV::loglik_calc()
{
	double loglik = 0.0;

  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += lgamma(beta + (double)n_kv(k, v) ) - lgamma(beta);
    }
    // word normalization
    loglik += lgamma( beta * (double)num_vocab ) - lgamma(beta * (double)num_vocab + (double)n_k(k) );
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

void LDACOV::lambda_store()
{
	Rcpp::NumericMatrix Lambda_R = Rcpp::wrap(Lambda);
	List Lambda_iter = model["Lambda_iter"];
	Lambda_iter.push_back(Lambda_R);
	model["Lambda_iter"] = Lambda_iter;
}

int LDACOV::sample_z(Eigen::VectorXd &alpha, int &z,
				  	     int &w, int &doc_id)
{
	// Remove data
	n_kv(z, w) -= 1;
	n_k(z) -= 1;
	n_dk(doc_id, z) -= 1;

	VectorXd z_prob_vec = VectorXd::Zero(num_topics);
	int new_z = -1;
	double numerator, denominator;

	for(int k=0; k<num_topics; k++){
		numerator = (beta + (double)n_kv(k, w)) *
			((double)n_k(k) + gamma_2) *
			((double)n_dk(doc_id, k) + alpha(k));

		denominator = ((double)num_vocab * beta + (double)n_k(k)) *
			((double)n_k(k) + gamma_1 + (double)n_k(k) + gamma_2);

		z_prob_vec(k) = numerator / denominator;	
	}

	double sum = z_prob_vec.sum();
	new_z = rcat_without_normalize(z_prob_vec, sum); // take a sample

  // add back data counts
	n_kv(new_z, w) += 1;
	n_k(new_z) += 1;
  n_dk(doc_id, new_z) += 1;

	return new_z;
}

int LDACOV::rcat_without_normalize(Eigen::VectorXd &prob,
		double &total)
{ // was called 'multi1'
  double u = R::runif(0, 1) * total;
  double temp = 0.0;
  int index = 0;
  for (int ii = 0; ii < prob.size(); ii++){
    temp += prob(ii);
    if (u < temp){
      index = ii;
      break;
    }
  }
  return index;
}

vector<int> LDACOV::shuffled_indexes(int m) {
  vector<int> v(m);
  iota(v.begin(), v.end(), 0);
  random_shuffle(v.begin(), v.end(), sampler::rand_wrapper);
  return v;
}


double LDACOV::expand(double& p, double& A){
	double res = -(1.0/A) * log((1.0/p) - 1.0);
	return res;
}

double LDACOV::shrink(double& x, double& A){
	double res = 1.0 / (1.0 + exp(-A*x));
	return res;
}
