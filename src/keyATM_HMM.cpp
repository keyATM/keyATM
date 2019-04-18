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

	// Initialize variables we use in the sampling
	logfy = VectorXd::Zero(num_states);
	st_k = VectorXd::Zero(num_states);
	logst_k = VectorXd::Zero(num_states);
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
	sample_forward();  // calculate Psk
	sample_backward();  // sample S_est
	sample_P();  // sample P_est

	store_S_est();
}

void keyATMhmm::sample_forward()
{ // Calculate Psk (num_doc, num_states)
	for(int d=0; d<num_doc; d++){
		if(d == 0){
			// First document should be the first state
			Psk(0, 0) = 1.0;
			continue;
		}	

		// Prepare f in Eq.(6) of Chib (1998)
		for(int s=0; s<num_states; s++){
			alpha = alphas.row(s).transpose();
			logfy(s) = polyapdfln(d, alpha);
		}	

		// Prepare Pst
		st_1l = Psk.row(d-1);  // previous observation
		st_k = (st_1l.transpose() * P_est);

		// Format numerator and calculate denominator at the same time
		logsum = 0.0;
		for(int s=0; s<num_states; s++){
			if(st_k(s) != 0.0){
				loglik = log(st_k(s)) + logfy(s);
				logst_k(s) = log(st_k(s)) + logfy(s);
				logsum = logsumexp(logsum, loglik, (s == 0));
			}
		}

		for(int s=0; s<num_states; s++){
			if(st_k(s) != 0.0){
				Psk(d, s) = exp( logst_k(s) - logsum );	
			}else{
				Psk(d, s) = 0.0;	
			}	
		}

	}

}

double keyATMhmm::polyapdfln(int &d, VectorXd &alpha)
{ // Polya distribution log-likelihood
	loglik = 0.0;

	loglik += lgamma( alpha.sum() ) - lgamma( doc_each_len[d] + alpha.sum() );
	for (int k = 0; k < num_topics; k++){
		loglik += lgamma( n_dk(d,k) + alpha(k) ) - lgamma( alpha(k) );
	}

	return loglik;
}


void keyATMhmm::sample_backward()
{
}


void keyATMhmm::sample_P()
{
}


void keyATMhmm::store_S_est()
{
	// Store state
	Rcpp::NumericVector state_R = Rcpp::wrap(S_est);
	List S_iter = model["S_iter"];
	S_iter.push_back(state_R);
	model["S_iter"] = S_iter;
}


double keyATMhmm::loglik_total()
{
  loglik = 0.0;
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
		
    loglik += lgamma( alpha.sum() ) - lgamma( doc_each_len[d] + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += lgamma( n_dk(d,k) + alpha(k) ) - lgamma( alpha(k) );
    }
  }

	return loglik;
}

