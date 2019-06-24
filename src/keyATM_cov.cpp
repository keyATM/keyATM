#include "keyATM_cov.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

keyATMcov::keyATMcov(List model_, const int iter_, const int output_per_) :
	keyATMbase(model_, iter_, output_per_) // pass to parent!
{
	// Constructor
	read_data();
	initialize();
	iteration();
	
	// Additional function
	
	// Add Sampling Info
	Rcpp::IntegerVector sampling_info = Rcpp::wrap(mh_info);
	List sampling_info_list = model["sampling_info"];
	sampling_info_list.push_back(sampling_info);
	model["sampling_info"] = sampling_info_list;
}


void keyATMcov::read_data_specific()
{
	// Covariate
	NumericMatrix C_r = model["C"];
	C = Rcpp::as<Eigen::MatrixXd>(C_r);
	num_cov = C.cols();

}

void keyATMcov::initialize_specific()
{
	// Alpha
	Alpha = MatrixXd::Zero(num_doc, num_topics);
	alpha = VectorXd::Zero(num_topics); // use in iteration

	// Lambda
	Lambda = MatrixXd::Zero(num_topics, num_cov);
	for(int k=0; k<num_topics; k++){
		// Initialize with R random
		for(int i=0; i<num_cov; i++){
			Lambda(k, i) = R::rnorm(0.0, 0.3);
		}
	}

}



void keyATMcov::iteration_single(int &it)
{ // Single iteration

	doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

	// Create Alpha for this iteration
	Alpha = (C * Lambda.transpose()).array().exp();

	for (int ii = 0; ii < num_doc; ii++){
		doc_id_ = doc_indexes[ii];
		doc_x = X[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
		doc_length = doc_each_len[doc_id_];
		
		token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
		
		// Prepare Alpha for the doc
		alpha = Alpha.row(doc_id_).transpose(); // take out alpha
		
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


void keyATMcov::sample_parameters()
{
	sample_lambda();
}


void keyATMcov::sample_lambda()
{
	// Sampling
	u = unif_rand(); // select sampling methods randomly

	// if(u < 0.3){
	// 	sample_lambda_mh();
	// }else if (u < 0.6){
	// 	sample_lambda_mh_single();
	// }else{
	// 	sample_lambda_slice();	
	// }

	if(u < 0.4){
		sample_lambda_mh_single();
	}else{
		sample_lambda_slice();
	}

	// Store Lambda
	Rcpp::NumericMatrix Lambda_R = Rcpp::wrap(Lambda);
	List Lambda_iter = model["Lambda_iter"];
	Lambda_iter.push_back(Lambda_R);
	model["Lambda_iter"] = Lambda_iter;
}


void keyATMcov::sample_lambda_mh_single()
{
	topic_ids = sampler::shuffled_indexes(num_topics);
	cov_ids = sampler::shuffled_indexes(num_cov);
	Lambda_current = 0.0;
	llk_current = 0.0;
	llk_proposal = 0.0;
	diffllk = 0.0;
	r = 0.0; u = 0.0;
	int k, t;

	for(int kk=0; kk<num_topics; kk++){
		k = topic_ids[kk];
		
		for(int tt=0; tt<num_cov; tt++){
			t = cov_ids[tt];
		
			mh_info[1] += 1; // how many times we run mh

			Lambda_current = Lambda(k,t);
			
			// Current llk
			llk_current = likelihood_lambda();
			
			// Proposal
			Lambda(k, t) += R::rnorm(0.0, mh_sigma);
			llk_proposal = likelihood_lambda();
			
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


void keyATMcov::sample_lambda_mh()
{
	
	topic_ids = sampler::shuffled_indexes(num_topics);
	VectorXd Lambda_current;
	llk_current = 0.0;
	llk_proposal = 0.0;
	diffllk = 0.0;
	r = 0.0; u =0.0;
	int k;

	for(int kk=0; kk<num_topics; kk++){
		k = topic_ids[kk];
		mh_info[1] += 1; // how many times we run mh

		Lambda_current = Lambda.row(k).transpose();
		
		// Current llk
		llk_current = likelihood_lambda();
		
		// Proposal
		proposal_lambda(k);
		llk_proposal = likelihood_lambda();
		
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


double keyATMcov::likelihood_lambda()
{
	double loglik = 0.0;
	Alpha = (C * Lambda.transpose()).array().exp();
	alpha = VectorXd::Zero(num_topics);

	for(int d=0; d<num_doc; d++){
		alpha = Alpha.row(d).transpose(); // Doc alpha, column vector
		// alpha = ((C.row(d) * Lambda)).array().exp(); // Doc alpha, column vector
	
		loglik += mylgamma(alpha.sum()); 
				// the first term numerator in the first square bracket
		loglik -= mylgamma( doc_each_len[d] + alpha.sum() ); 
				// the second term denoinator in the first square bracket
	
		for(int k=0; k<num_topics; k++){
			loglik -= mylgamma(alpha(k));
				// the first term denominator in the first square bracket
			loglik += mylgamma( n_dk(d, k) + alpha(k) );
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



void keyATMcov::proposal_lambda(int& k){

	for(int i=0; i<num_cov; i++){
		Lambda.coeffRef(k, i) += R::rnorm(0.0, mh_sigma);
	}

}


void keyATMcov::sample_lambda_slice()
{
	// Slice sampling for Lambda

	start = 0.0; end = 0.0;
	previous_p = 0.0; new_p = 0.0;
	newlikelihood = 0.0; slice_ = 0.0; current_lambda = 0.0;
  topic_ids = sampler::shuffled_indexes(num_topics);
	cov_ids = sampler::shuffled_indexes(num_cov);
	int k, t;
	const double A = slice_A;

	store_loglik = likelihood_lambda();
	newlambdallk = 0.0;

	for(int kk=0; kk<num_topics; kk++){
		k = topic_ids[kk];

		for(int tt=0; tt<num_cov; tt++){
			t = cov_ids[tt];

			start = 0.0; // shrink
			end = 1.0; // shrink

			current_lambda = Lambda(k,t);
			previous_p = shrink(current_lambda, A);
			slice_ = store_loglik - std::log(A * previous_p * (1.0 - previous_p)) 
							+ log(unif_rand()); // <-- using R random uniform


			for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++){
				new_p = sampler::slice_uniform(start, end); // <-- using R function above
				Lambda(k,t) = expand(new_p, A); // expand

				newlambdallk = likelihood_lambda();

				newlikelihood = newlambdallk - std::log(A * new_p * (1.0 - new_p));

				if (slice_ < newlikelihood){
					store_loglik = newlambdallk;
					break;
				} else if (previous_p < new_p){
					end = new_p;
				} else if (new_p < previous_p){
					start = new_p;
				} else {
					// Rcerr << "Something goes wrong in sample_lambda_slice()" << std::endl;
					Rcpp::stop("Something goes wrong in sample_lambda_slice(). Adjust `A_slice`.");
					Lambda(k,t) = current_lambda;
					break;
				}

			} // for loop for shrink time

		} // for loop for num_cov
	} // for loop for num_topics
}


double keyATMcov::loglik_total()
{
  double loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += mylgamma(beta + n_x0_kv(k, v) / vocab_weights(v) ) - mylgamma(beta);
      // loglik += mylgamma(beta_s + n_x1_kv.coeffRef(k, v) / vocab_weights(v) ) - mylgamma(beta_s);
    }

		// n_x1_kv
		for (SparseMatrix<double,RowMajor>::InnerIterator it(n_x1_kv, k); it; ++it){
			loglik += mylgamma(beta_s + it.value() / vocab_weights(it.index()) ) - mylgamma(beta_s);
		}

    // word normalization
    loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_x0_k_noWeight(k) );
    loglik += mylgamma( beta_s * (double)num_vocab ) - mylgamma(beta_s * (double)num_vocab + n_x1_k_noWeight(k) );
    // x
    loglik += mylgamma( n_x0_k_noWeight(k) + gamma_2 ) - mylgamma(n_x1_k_noWeight(k) + gamma_1 + n_x0_k_noWeight(k) + gamma_2)
      + mylgamma( n_x1_k_noWeight(k) + gamma_1 ) ;
		
    // x normalization
    loglik += mylgamma(gamma_1 + gamma_2) - mylgamma(gamma_1) - mylgamma(gamma_2);
  }


  // z
	Alpha = (C * Lambda.transpose()).array().exp();
	alpha = VectorXd::Zero(num_topics);

  for (int d = 0; d < num_doc; d++){
		alpha = Alpha.row(d).transpose(); // Doc alpha, column vector	
		
    loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len[d] + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }
  }

  return loglik;
}
