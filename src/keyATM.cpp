#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

keyATMbase::keyATMbase(List model_, const int iter_, const int output_per_)
{
	model = model_;
	iter = iter_;
	output_per = output_per_;
}

void keyATMbase::read_data()
{
	read_data_common();
	read_data_specific();
}

void keyATMbase::read_data_common()
{
	// Read data
	W = model["W"], Z = model["Z"], X = model["X"];
  files = model["files"], vocab = model["vocab"];
  nv_alpha = model["alpha"];
  gamma_1 = model["gamma_1"], gamma_2 = model["gamma_2"];
  beta = model["beta"], beta_s = model["beta_s"];
  k_free = model["extra_k"];
  seeds = model["seeds"];
  k_seeded = seeds.size();
	model_fit = model["model_fit"];

  num_topics = k_seeded + k_free;
  // alpha -> specific function	

}


void keyATMbase::initialize()
{
	initialize_common();
	initialize_specific();
}


void keyATMbase::initialize_common()
{
	// Vector that stores seed words (words in dictionary)
  for (int ii = 0; ii < k_seeded; ii++){
    IntegerVector wd_ids = seeds[ii];
    seed_num.push_back(wd_ids.size());
		
    std::unordered_set<int> keywords_set;
    for (int jj = 0; jj < wd_ids.size(); jj++){
			int wd_id = wd_ids(jj);
			keywords_set.insert(wd_id);
		}
		
    keywords.push_back(keywords_set);
  }

	for(int i=k_seeded; i<num_topics; i++){
		std::unordered_set<int> keywords_set{ -1 };
	
		seed_num.push_back(0);
		keywords.push_back(keywords_set);
	}

  // document-related constants
  num_vocab = vocab.size();
	num_doc = files.size();

  // storage for sufficient statistics and their margins
  n_x0_kv = MatrixXi::Zero(num_topics, num_vocab);
  n_x1_kv = MatrixXi::Zero(num_topics, num_vocab);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_x0_k = VectorXi::Zero(num_topics);
  n_x1_k = VectorXi::Zero(num_topics);

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
  total_words = (int)n_dk.sum();
	
	prepare_data = std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::high_resolution_clock::now() - start).count();
}


// Iteration
void keyATMbase::iteration()
{
	for(int it=0; it<iter; it++){
		iteration_single();	

		int r_index = it + 1;
		if(r_index % output_per == 0 || r_index == 1 || r_index == iter ){
			sampling_store(r_index);
		}

		checkUserInterrupt();
	}

	model["model_fit"] = model_fit;
}

void keyATMbase::sampling_store(int &r_index)
{
	double loglik = loglik_total();
	double perplexity = exp(-loglik / (double)total_words);

	NumericVector model_fit_vec;
	model_fit_vec.push_back(r_index);
	model_fit_vec.push_back(loglik);
	model_fit_vec.push_back(perplexity);
	model_fit.push_back(model_fit_vec);

	Rcerr << "[" << r_index << "] log likelihood: " << loglik <<
					 " (perplexity: " << perplexity << ")" << std::endl;
}

// Sampling
int keyATMbase::sample_z(VectorXd &alpha, int &z, int &x,
												 int &w, int &doc_id)
{
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


int keyATMbase::sample_x(VectorXd &alpha, int &z, int &x,
				  	     int &w, int &doc_id)
{
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



// Utilities
double keyATMbase::logsumexp(double &x, double &y, bool flg){
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


double keyATMbase::logsumexp_Eigen(VectorXd &vec){
  double sum = 0.0;
  for(int i = 0; i < vec.size(); i++){
    sum = logsumexp(sum, vec[i], (i == 0));
  }
  return sum;
}


double keyATMbase::gammapdfln(const double x, const double a, const double b){
  return a * log(b) - lgamma(a) + (a - 1.0) * log(x) - b * x;
}


NumericVector keyATMbase::alpha_reformat(VectorXd& alpha, int& num_topics){
	NumericVector alpha_rvec(num_topics);

	for(int i=0; i<num_topics; ++i){
		alpha_rvec[i] = alpha(i);
	}

	return alpha_rvec;
}


double keyATMbase::expand(double& p){
	double res = -(1.0/slice_A) * log((1.0/p) - 1.0);
	return res;
}

double keyATMbase::shrink(double& x){
	double res = 1.0 / (1.0 + exp(-slice_A*x));
	return res;
}


List keyATMbase::return_model(){
	return model;
}
