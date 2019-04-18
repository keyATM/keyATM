#include "keyATM_HMM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

keyATMhmm::keyATMhmm(List model_, const int iter_, const int output_per_) :
  keyATMbase(model_, iter_, output_per_) // pass to parent!
{
	// Constructor
	read_data();
	initialize();
	iteration();
}


void keyATMhmm::read_data_specific()
{
	num_states = model["num_states"];
	index_states = num_states - 1;
}


void keyATMhmm::initialize_specific()
{
	// Initialize Psk
	Psk = MatrixXd::Zero(num_doc, num_states);

	// Initialize S_est
	// Use multinomial distribution (with flat probability)
	// to decide the number of each state
	// and push it into S_est.
	VectorXi S_est_num = VectorXi::Zero(num_states);
	VectorXd S_est_temp = VectorXd::Zero(num_states);
	double cumulative = 1.0 / num_states;
	double u;
	int index;
	for(int i=0; i<num_states; i++){
		S_est_temp(i) = cumulative * (i+1);
	}

	for(int j=0; j<num_doc; j++){
		u = R::runif(0, 1);
		for(int i=0; i<num_states; i++){
			if(u < S_est_temp(i)){
				index = i;
				break;
			}
		}
		S_est_num(index) += 1;
	}

	S_est = VectorXi::Zero(num_doc);
	int count;
	index = 0;
	for(int i=0; i<num_states; i++){
		count = S_est_num(i);
		for(int j=0; j<count; j++){
			S_est(index) = i;
			index += 1;
		}
	}

	// Initializae P_est
	P_est = MatrixXd::Zero(num_states, num_states);
	double prob;
	for(int i=0; i<=(index_states-1); i++){
		prob = R::rbeta(1.0, 1.0);
		P_est(i, i) = prob;
		P_est(i, i+1) = 1-prob;
	}
	P_est(index_states, index_states) = 1;


	// Initialize alphas;
	alphas = MatrixXd::Constant(num_states, num_topics, 50.0/num_topics);
}


void keyATMhmm::iteration_single()
{
	doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

	for (int ii = 0; ii < num_doc; ii++){
		doc_id_ = doc_indexes[ii];
		doc_x = X[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
		doc_length = doc_each_len[doc_id_];

		alpha = alphas.row(S_est(doc_id_)).transpose(); // select alpha for this document
		
		token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
		
		// Iterate each word in the document
		for (int jj = 0; jj < doc_length; jj++){
			w_position = token_indexes[jj];
			x_ = doc_x[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position];
		
			new_z = sample_z(alpha, z_, x_, w_, doc_id_);
			doc_z[w_position] = new_z;
		
	
			z_ = doc_z[w_position]; // use updated z
			new_x = sample_x(alpha, z_, x_, w_, doc_id_);
			doc_x[w_position] = new_x;
		}
		
		Z[doc_id_] = doc_z;
		X[doc_id_] = doc_x;
	}
	sample_parameters();
}


void keyATMhmm::sample_parameters()
{

}


double keyATMhmm::loglik_total()
{
  double loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += lgamma(beta + n_x0_kv(k, v) ) - lgamma(beta);
      loglik += lgamma(beta_s + n_x1_kv(k, v) ) - lgamma(beta_s);
    }
    // word normalization
    loglik += lgamma( beta * (double)num_vocab ) - lgamma(beta * (double)num_vocab + n_x0_k(k) );
    loglik += lgamma( beta_s * (double)num_vocab ) - lgamma(beta_s * (double)num_vocab + n_x1_k(k) );
    // x
    loglik += lgamma( n_x0_k(k) + gamma_2 ) - lgamma(n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2)
      + lgamma( n_x1_k(k) + gamma_1 ) ;
		
    // x normalization
    loglik += lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) - lgamma(gamma_2);
  }


  // z
  for (int d = 0; d < num_doc; d++){
		alpha = alphas.row(S_est(doc_id_)).transpose(); // Doc alpha, column vector	
		
    loglik += lgamma( alpha.sum() ) - lgamma( n_dk.row(d).sum() + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += lgamma( n_dk(d,k) + alpha(k) ) - lgamma( alpha(k) );
    }
  }
}

