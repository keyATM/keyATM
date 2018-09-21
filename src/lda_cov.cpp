#include "lda_cov.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

namespace sampler{
	inline int rand_wrapper(const int n) { return floor(unif_rand() * n); }
}

LDACOV::LDACOV(List model, const int K, const int iter_, const int output_iter_)
{
	// Get Info from R
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
	Eigen::VectorXd alpha_ = VectorXd::Zero(num_topics);
	
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
	
	}
}

int LDACOV::sample_z(Eigen::VectorXd &alpha, int &z,
				  	     int &w, int &doc_id)
{
	return 1;
}

vector<int> LDACOV::shuffled_indexes(int m) {
  vector<int> v(m);
  iota(v.begin(), v.end(), 0);
  random_shuffle(v.begin(), v.end(), sampler::rand_wrapper);
  return v;
}

