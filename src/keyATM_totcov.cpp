#include "keyATM_totcov.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

keyATMtotcov::keyATMtotcov(List model_, const int iter_, const int output_per_) :
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


void keyATMtotcov::read_data_specific()
{
	// Covariate
	NumericMatrix C_r = model["C"];
	C = Rcpp::as<Eigen::MatrixXd>(C_r);
	num_cov = C.cols();

	// TOT
	// Read time stamp
	timestamps = Rcpp::as<Eigen::VectorXd>(model["timestamps"]);

	// Options
	List options = model["options"];
	logsumexp_approx = options["logsumexp_approx"];
	use_mom = options["use_mom"];
}


void keyATMtotcov::initialize_specific()
{
	// Covariate
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


	// TOT
	// Store parameters for time Beta
	beta_params = MatrixXd::Constant(num_topics, 2, 0.5);

	beta_tg = VectorXd::Zero(num_topics);
	beta_lg = VectorXd::Zero(num_topics);
	beta_tg_base = VectorXd::Zero(num_topics);
	beta_lg_base = VectorXd::Zero(num_topics);

	// Store time stamps
	for(int k=0; k<num_topics; k++){
		vector<double> temp;
		store_t.push_back(temp);
	}
}


void keyATMtotcov::iteration_single(int &it)
{// Single iteration

	doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

	// Create Alpha for this iteration
	Alpha = (C * Lambda.transpose()).array().exp();

	// Clear temporary time stamp
	for(int k=0; k<num_topics; k++){
		store_t[k].clear();
	}

	// Apply gamma
	for(int k=0; k<num_topics; k++){ 
		beta_a = beta_params(k, 0);
		beta_b = beta_params(k, 1);	

		// Log version
		if(!logsumexp_approx)
			beta_lg_base(k) =  mylgamma(beta_a + beta_b) - (mylgamma(beta_a) + mylgamma(beta_b));	

		// Normal version
		beta_tg_base(k) =  tgamma(beta_a + beta_b) / (tgamma(beta_a) * tgamma(beta_b));
	}


	// Sampling
	for (int ii = 0; ii < num_doc; ii++){
		doc_id_ = doc_indexes[ii];
		doc_x = X[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
		doc_length = doc_each_len[doc_id_];
		
		token_indexes = sampler::shuffled_indexes(doc_length); //shuffle

		// Prepare beta_a and beta_b for sampling
		timestamp_d = timestamps(doc_id_);

		use_log = 0;
		for(int k=0; k<num_topics; k++){ 
			beta_a = beta_params(k, 0);
			beta_b = beta_params(k, 1);
		
			// for log version
			if(!logsumexp_approx)
				beta_lg(k) = beta_lg_base(k) + (beta_a - 1.0) * log(1.0 - timestamp_d) +
										 (beta_b - 1.0) * log(timestamp_d);
		
			// For normal version
			check_frac = beta_tg_base(k) * pow(1.0 - timestamp_d, beta_a - 1.0) * pow(timestamp_d, beta_b - 1.0);

			if(check_frac < numeric_limits<double>::min() || 
					check_frac > numeric_limits<double>::max()){
				if(logsumexp_approx){
					beta_tg(k) = 2.22507e-200;	
				}else{
					use_log = 1;
				}
			}else{
				beta_tg(k) = check_frac;
			}

		}

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



int keyATMtotcov::sample_z_log(VectorXd &alpha, int &z, int &x,
																		 int &w, int &doc_id)
{
  // removing data is already done in sample_z

  new_z = -1; // debug
  if (x == 0){
    for (int k = 0; k < num_topics; ++k){

      numerator = log( (beta + n_x0_kv(k, w)) *
        (n_x0_k(k) + gamma_2) *
        (n_dk(doc_id, k) + alpha(k)) );

      denominator = log( ((double)num_vocab * beta + n_x0_k(k)) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2) );

      z_prob_vec(k) = numerator - denominator + beta_lg(k);
    }


		sum = logsumexp_Eigen(z_prob_vec, num_topics);
		z_prob_vec = (z_prob_vec.array() - sum).exp();
		new_z = sampler::rcat(z_prob_vec, num_topics);


  } else {
    for (int k = 0; k < num_topics; ++k){
      if (keywords[k].find(w) == keywords[k].end()){
        z_prob_vec(k) = 0.0;
				continue;
      } else{ 

        numerator = log( (beta_s + n_x1_kv.coeffRef(k, w)) *
           (n_x1_k(k) + gamma_1)  *
           (n_dk(doc_id, k) + alpha(k)) );
      }
      denominator = log( ((double)seed_num[k] * beta_s + n_x1_k(k) ) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2) );

      z_prob_vec(k) = numerator - denominator + beta_lg(k);
    }



		sum = logsumexp_Eigen(z_prob_vec, num_topics);
		z_prob_vec = (z_prob_vec.array() - sum).exp();
		new_z = sampler::rcat(z_prob_vec, num_topics);

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

	// Store time stamp
	store_t[new_z].push_back(timestamp_d);

  return new_z;
}


int keyATMtotcov::sample_z(VectorXd &alpha, int &z, int &x,
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

	if(use_log == 1)
		return sample_z_log(alpha, z, x, w, doc_id);

  if (x == 0){
    for (int k = 0; k < num_topics; ++k){
      numerator = (beta + n_x0_kv(k, w)) *
        (n_x0_k(k) + gamma_2) *
        (n_dk(doc_id, k) + alpha(k));

      denominator = ((double)num_vocab * beta + n_x0_k(k)) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2);


			check_frac = numerator / denominator * beta_tg(k);

			if(check_frac < numeric_limits<double>::min() || 
					check_frac > numeric_limits<double>::max()){
					cout << check_frac << " ";  // debug
				if(logsumexp_approx){
					check_frac = 2.22507e-200;	
				}else{
					return sample_z_log(alpha, z, x, w, doc_id);
				}
			}

			z_prob_vec(k) = check_frac;
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
           (n_x1_k(k) + gamma_1)  *
           (n_dk(doc_id, k) + alpha(k));
      }
      denominator = ((double)seed_num[k] * beta_s + n_x1_k(k) ) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2);


			check_frac = numerator / denominator * beta_tg(k);

			if(check_frac < numeric_limits<double>::min() | 
					check_frac > numeric_limits<double>::max()){
				return sample_z_log(alpha, z, x, w, doc_id);
			}

			z_prob_vec(k) = check_frac;
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

	// Store time stamp
	store_t[new_z].push_back(timestamp_d);

  return new_z;
}


void keyATMtotcov::sample_parameters()
{
	sample_lambda();
	sample_betaparam();
}




void keyATMtotcov::sample_lambda()
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


void keyATMtotcov::sample_lambda_mh_single()
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


void keyATMtotcov::sample_lambda_mh()
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


double keyATMtotcov::likelihood_lambda()
{
	double loglik = 0.0;
	Alpha = (C * Lambda.transpose()).array().exp();
	alpha = VectorXd::Zero(num_topics);

	for(int d=0; d<num_doc; d++){
		alpha = Alpha.row(d).transpose(); // Doc alpha, column vector
	
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



void keyATMtotcov::proposal_lambda(int& k){

	for(int i=0; i<num_cov; i++){
		Lambda.coeffRef(k, i) += R::rnorm(0.0, mh_sigma);
	}

}


void keyATMtotcov::sample_lambda_slice()
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



void keyATMtotcov::sample_betaparam(){

	if(use_mom){
		for(int k=0; k<num_topics; k++){
			if(store_t[k].size() == 0)
				continue;
		
			timestamps_k = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
														(store_t[k].data(), store_t[k].size());		
		
			beta_mean = timestamps_k.mean();
			beta_var = ( (timestamps_k.array() - beta_mean) * (timestamps_k.array() - beta_mean) ).sum() /
										(store_t[k].size() - 1.0);

			if(beta_var == 0.0)
				continue;
		
			beta_var = 1.0 / (beta_var);  // beta_var reciprocal
			beta_var = ( (beta_mean * (1-beta_mean)) * beta_var - 1.0);
			
			beta_params(k, 0) = beta_mean * beta_var;
			beta_params(k, 1) = (1.0 - beta_mean) * beta_var;
		}

	}else{

		topic_ids = sampler::shuffled_indexes(num_topics);

		for(int j=0; j<num_topics; j++){
			for(int i=0; i<2; i++){
				k = topic_ids[j];
				current_param = beta_params(k, i);

				// start = min_v / (1.0 + min_v);
				// end = 1.0;

				start = min_v / (1.0 + min_v);
				end = 15.0 / (1.0 + 15.0);

				previous_p = current_param / (1.0 + current_param);

				temp_beta_loglik = beta_loglik(k, i);

				if(temp_beta_loglik == -1.0)
					continue;

				slice_ = temp_beta_loglik - 2.0 * log(1.0 - previous_p) 
									+ log(unif_rand()); 


				for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++){
					new_p = sampler::slice_uniform(start, end); // <-- using R function above
					beta_params(k, i) = new_p / (1.0 - new_p); // expandp

					newlikelihood = beta_loglik(k, i) - 2.0 * log(1.0 - new_p);

					if (slice_ < newlikelihood){
						break;
					} else if (previous_p < new_p){
						end = new_p;
					} else if (new_p < previous_p){
						start = new_p;
					} else {
						Rcpp::stop("Something goes wrong in sample_lambda_slice(). Adjust `A_slice`.");
						beta_params(k, i) = current_param;
						break;
					}
				}
			
			}  // for i	
		}  // for j

	}
}


double keyATMtotcov::beta_loglik(const int &k, const int &i){
	loglik = 0.0;

  for (int d = 0; d < num_doc; d++){
		loglik += n_dk(d, k) * betapdfln(timestamps(d), beta_params(k,0), beta_params(k,1));
  }

	if(loglik == 0)
		return -1.0;

	// Prior
	loglik += gammapdfln(beta_params(k, i), ts_g1, ts_g2);

	return loglik;

}


double keyATMtotcov::loglik_total()
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

			// time stamps
			loglik += n_dk(d, k) * betapdfln(timestamps(d), beta_params(k,0), beta_params(k,1));
    }
  }

	// Time stamps Prior
	for(int k=0; k<num_topics; k++){
		for(int i=0; i<2; i++){
			loglik += gammapdfln(beta_params(k, i), ts_g1, ts_g2);
		}	
	}

  return loglik;

}






