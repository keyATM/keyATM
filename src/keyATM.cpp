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

keyATMbase::~keyATMbase(){

}

void keyATMbase::read_data()
{
	read_data_common();
	read_data_specific();
}

void keyATMbase::read_data_common()
{
	// Read data
	W = model["W"]; Z = model["Z"]; X = model["X"];
  files = model["files"]; vocab = model["vocab"];
  nv_alpha = model["alpha"];
  gamma_1 = model["gamma_1"];
	gamma_2 = model["gamma_2"];
  beta = model["beta"];
	beta_s = model["beta_s"];
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

	// prepare_data = std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::high_resolution_clock::now() - start).count();
}


void keyATMbase::initialize_common()
{
	// Vector that stores seed words (words in dictionary)
	int wd_id;
	IntegerVector wd_ids;
  for (int ii = 0; ii < k_seeded; ii++){
    wd_ids = seeds[ii];
    seed_num.push_back(wd_ids.size());
		
    std::unordered_set<int> keywords_set;
    for (int jj = 0; jj < wd_ids.size(); jj++){
			wd_id = wd_ids(jj);
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
  // n_x0_kv.resize(num_topics, num_vocab);
	n_x0_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_x1_kv.resize(num_topics, num_vocab);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_x0_k = VectorXd::Zero(num_topics);
  n_x0_k_noWeight = VectorXd::Zero(num_topics);
  n_x1_k = VectorXd::Zero(num_topics);
  n_x1_k_noWeight = VectorXd::Zero(num_topics);
	vocab_weights = VectorXd::Constant(num_vocab, 1.0);

	int x, z, w;
	int doc_len;
	IntegerVector doc_x, doc_z, doc_w;


	// Construct vocab weights
  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    doc_w = W[doc_id];
		doc_len = doc_w.size();
		doc_each_len.push_back(doc_len);
	
    for(int w_position = 0; w_position < doc_len; w_position++){
			w = doc_w[w_position];
			vocab_weights(w) += 1.0;
    }
  }
  total_words = (int)vocab_weights.sum();
	vocab_weights = vocab_weights.array() / (double)total_words;
	vocab_weights = vocab_weights.array().log();
	vocab_weights = - vocab_weights.array() / log(2);
	

	if(use_weight == 0){
		cout << "Not using weights!! Check use_weight" << endl;
		vocab_weights = VectorXd::Constant(num_vocab, 1.0);
	}


	// Construct data matrices
	vector<Triplet> trip_x1; 	
	
  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
		doc_len = doc_each_len[doc_id];

    for(int w_position = 0; w_position < doc_len; w_position++){
      x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];
      if (x == 0){
        n_x0_kv(z, w) += vocab_weights(w);
				// trip_x0.push_back(Triplet(z, w, vocab_weights(w)));
        n_x0_k(z) += vocab_weights(w);
        n_x0_k_noWeight(z) += 1.0;
      } else {
        // n_x1_kv(z, w) += vocab_weights(w);
				trip_x1.push_back(Triplet(z, w, vocab_weights(w)));
        n_x1_k(z) += vocab_weights(w);
        n_x1_k_noWeight(z) += 1.0;
      }
      n_dk(doc_id, z) += 1.0;
    }
  }
	// n_x0_kv.setFromTriplets(trip_x0.begin(), trip_x0.end());
	n_x1_kv.setFromTriplets(trip_x1.begin(), trip_x1.end());
	

	// Use during the iteration
	z_prob_vec = VectorXd::Zero(num_topics);
	
}


// Iteration
void keyATMbase::iteration()
{
	for(int it=0; it<iter; it++){
		iteration_single(it);	

		int r_index = it + 1;
		if(r_index % output_per == 0 || r_index == 1 || r_index == iter ){
			sampling_store(r_index);
			verbose_special(r_index);
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

void keyATMbase::verbose_special(int &r_index){
	// If there is anything special to show or store, write here.
}

// Sampling
int keyATMbase::sample_z(VectorXd &alpha, int &z, int &x,
												 int &w, int &doc_id)
{
  // remove data
  if (x == 0){
    n_x0_kv(z, w) -= vocab_weights(w);
    n_x0_k(z) -= vocab_weights(w);
    n_x0_k_noWeight(z) -= 1.0;
  } else if (x==1) {
    n_x1_kv.coeffRef(z, w) -= vocab_weights(w);
    n_x1_k(z) -= vocab_weights(w);
    n_x1_k_noWeight(z) -= 1.0;
  } else {
		Rcerr << "Error at sample_z, remove" << std::endl;
	}

  n_dk(doc_id, z) -= 1;

  new_z = -1; // debug
  if (x == 0){
    for (int k = 0; k < num_topics; ++k){

      numerator = (beta + n_x0_kv(k, w)) *
        (n_x0_k(k) + gamma_2) *
        (n_dk(doc_id, k) + alpha(k));

      denominator = ((double)num_vocab * beta + n_x0_k(k)) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2);

      z_prob_vec(k) = numerator / denominator;
    }

    sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  } else {
    for (int k = 0; k < num_topics; ++k){
      if (keywords[k].find(w) == keywords[k].end()){
        z_prob_vec(k) = 0.0;
				continue;
      } else{ 
        numerator = (beta_s + n_x1_kv.coeffRef(k, w)) *
          (n_x1_k(k) + gamma_1) *
          (n_dk(doc_id, k) + alpha(k));
      }
      denominator = ((double)seed_num[k] * beta_s + n_x1_k(k) ) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2);

      z_prob_vec(k) = numerator / denominator;
    }


		sum = z_prob_vec.sum();
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  }

  // add back data counts
  if (x == 0){
    n_x0_kv(new_z, w) += vocab_weights(w);
    n_x0_k(new_z) += vocab_weights(w);
    n_x0_k_noWeight(new_z) += 1.0;
  } else if (x==1) {
    n_x1_kv.coeffRef(new_z, w) += vocab_weights(w);
    n_x1_k(new_z) += vocab_weights(w);
    n_x1_k_noWeight(new_z) += 1.0;
  } else {
		Rcerr << "Error at sample_z, add" << std::endl;
	}
  n_dk(doc_id, new_z) += 1;

  return new_z;
}


int keyATMbase::sample_x(VectorXd &alpha, int &z, int &x,
				  	     int &w, int &doc_id)
{
			
	// If a word is not a keyword, no need to sample
	if(keywords[z].find(w) == keywords[z].end())
		return x;
	
  // remove data
  if (x == 0){
    n_x0_kv(z, w) -= vocab_weights(w);
    n_x0_k(z) -= vocab_weights(w);
    n_x0_k_noWeight(z) -= 1.0;
  } else {
    n_x1_kv.coeffRef(z, w) -= vocab_weights(w);
    n_x1_k(z) -= vocab_weights(w);
    n_x1_k_noWeight(z) -= 1.0;
  }

  // newprob_x1()
  k = z;

	numerator = (beta_s + n_x1_kv.coeffRef(k, w)) *
		( n_x1_k(k) + gamma_1 );
	denominator = ((double)seed_num[k] * beta_s + n_x1_k(k) ) *
		(n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2);
	x1_prob = numerator / denominator;

	// newprob_x0()
	numerator = (beta + n_x0_kv(k, w)) *
		(n_x0_k(k) + gamma_2);

	denominator = ((double)num_vocab * beta + n_x0_k(k) ) *
		(n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2);
	x0_prob = numerator / denominator;

	// Normalize
	sum = x0_prob + x1_prob;

	x1_prob = x1_prob / sum;
	new_x = R::runif(0,1) <= x1_prob;  //new_x = Bern(x0_prob, x1_prob);

  // add back data counts
  if (new_x == 0){
    n_x0_kv(z, w) += vocab_weights(w);
    n_x0_k(z) += vocab_weights(w);
    n_x0_k_noWeight(z) += 1.0;
  } else {
    n_x1_kv.coeffRef(z, w) += vocab_weights(w);
    n_x1_k(z) += vocab_weights(w);
    n_x1_k_noWeight(z) += 1.0;
  }

  return new_x;
}



// Utilities
double keyATMbase::gammapdfln(const double x, const double a, const double b){
  return a * log(b) - lgamma(a) + (a - 1.0) * log(x) - b * x;
}


double keyATMbase::betapdf(const double x, const double a, const double b){
	return tgamma(a+b) / (tgamma(a) * tgamma(b)) * pow(x, a-1) * pow(1-x, b-1);
}

double keyATMbase::betapdfln(const double x, const double a, const double b){
	return (a-1)*log(x) + (b-1)*log(1.0-x) + mylgamma(a+b) - mylgamma(a) - mylgamma(b);
}

NumericVector keyATMbase::alpha_reformat(VectorXd& alpha, int& num_topics){
	NumericVector alpha_rvec(num_topics);

	for(int i=0; i<num_topics; ++i){
		alpha_rvec[i] = alpha(i);
	}

	return alpha_rvec;
}


double keyATMbase::gammaln_frac(const double &value, const int &count){
	// Calculate \log \frac{\gamma(value + count)}{\gamma(\value)}
	// Not very fast
	
	if(count > 19){
		return lgamma(value + count) - lgamma(value);	
	}else{
		gammaln_val = 0.0;

		for(int i=0; i<count; i++){
			gammaln_val += log(value + i);	
		}

		return gammaln_val;
	}
}


List keyATMbase::return_model(){
	return model;
}



