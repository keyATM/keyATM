#include "LDA_weightTOT.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


LDAweightTOT::LDAweightTOT(List model_, const int iter_, const int output_per_, const int use_weight_) :
	keyATMbase(model_, iter_, output_per_) // pass to parent!
{
	use_weight = use_weight_;
	// Constructor
	read_data();
	initialize();
	iteration();
}


void LDAweightTOT::read_data_specific()
{
	alpha = Rcpp::as<Eigen::VectorXd>(nv_alpha);

	// Read time stamp
	timestamps = Rcpp::as<Eigen::VectorXd>(model["timestamps"]);
}


void LDAweightTOT::initialize_specific()
{
	// Initialization for LDA weights 

  n_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_k = VectorXd::Zero(num_topics);
  n_k_noWeight = VectorXd::Zero(num_topics);

	int z, w;
	int doc_len;
	IntegerVector doc_z, doc_w;


	// Construct data matrices
  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    doc_z = Z[doc_id], doc_w = W[doc_id];
		doc_len = doc_each_len[doc_id];

    for(int w_position = 0; w_position < doc_len; w_position++){
      z = doc_z[w_position], w = doc_w[w_position];

			n_kv(z, w) += vocab_weights(w);
			n_k(z) += vocab_weights(w);
			n_k_noWeight(z) += 1.0;
      n_dk(doc_id, z) += 1.0;
    }
  }
	

	// Store parameters for time Beta
	beta_params = MatrixXd::Constant(num_topics, 2, 0.5);

	// Use during the iteration
	z_prob_vec = VectorXd::Zero(num_topics);

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

void LDAweightTOT::iteration_single(int &it)
{ // Single iteration


	// Clear temporary time stamp
	for(int k=0; k<num_topics; k++){
		store_t[k].clear();
	}

	// Apply gamma
	for(int k=0; k<num_topics; k++){ 
		beta_a = beta_params(k, 0);
		beta_b = beta_params(k, 1);	

		// Log version
		beta_lg_base(k) =  mylgamma(beta_a + beta_b) - (mylgamma(beta_a) + mylgamma(beta_b));	

		// Normal version
		beta_tg_base(k) =  tgamma(beta_a + beta_b) / (tgamma(beta_a) * tgamma(beta_b));
	}


	x_ = -1;  // we do not use x_ in LDA weight
	doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

	for (int ii = 0; ii < num_doc; ii++){
		doc_id_ = doc_indexes[ii];
		doc_z = Z[doc_id_], doc_w = W[doc_id_];
		doc_length = doc_each_len[doc_id_];
		
		token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
	
	
		// Prepare beta_a and beta_b for sampling
		timestamp_d = timestamps(doc_id_);
	
		use_log = 0;
		for(int k=0; k<num_topics; k++){ 
			beta_a = beta_params(k, 0);
			beta_b = beta_params(k, 1);
		
			// for log version
			beta_lg(k) = beta_lg_base(k) + (beta_a - 1.0) * log(1.0 - timestamp_d) + (beta_b - 1.0) * log(timestamp_d);
		
			// For normal version
			check_frac = beta_tg_base(k) * pow(1.0 - timestamp_d, beta_a - 1.0) * pow(timestamp_d, beta_b - 1.0);
			
			if(check_frac < numeric_limits<double>::min() | 
					check_frac > numeric_limits<double>::max()){
				use_log = 1;
			}else{
				beta_tg(k) = check_frac;
			}
			
		}
	
		
		// Iterate each word in the document
		for (int jj = 0; jj < doc_length; jj++){
			w_position = token_indexes[jj];
			z_ = doc_z[w_position], w_ = doc_w[w_position];
		
			new_z = sample_z(alpha, z_, x_, w_, doc_id_);
			doc_z[w_position] = new_z;
		}
		
		Z[doc_id_] = doc_z;
	}
	sample_parameters();

}


// Sampling
int LDAweightTOT::sample_z(VectorXd &alpha, int &z, int &x,
												 int &w, int &doc_id)
{
  // remove data
	n_kv(z, w) -= vocab_weights(w);
	n_k(z) -= vocab_weights(w);
	n_k_noWeight(z) -= 1.0;
  n_dk(doc_id, z) -= 1;

  new_z = -1; // debug

	if(use_log == 1)
		return sample_z_log(alpha, z, x, w, doc_id);

	for (int k = 0; k < num_topics; ++k){

		numerator = (beta + n_kv(k, w)) *
			(n_dk(doc_id, k) + alpha(k));

		denominator = ((double)num_vocab * beta + n_k(k)) ;

		check_frac = numerator / denominator * beta_tg(k);

		if(check_frac < numeric_limits<double>::min() | 
				check_frac > numeric_limits<double>::max()){
			return sample_z_log(alpha, z, x, w, doc_id);
		}

		z_prob_vec(k) = check_frac;
	}

	sum = z_prob_vec.sum(); // normalize
	new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample


  // add back data counts
	n_kv(new_z, w) += vocab_weights(w);
	n_k(new_z) += vocab_weights(w);
	n_k_noWeight(new_z) += 1.0;
  n_dk(doc_id, new_z) += 1;

	// Store time stamp
	store_t[new_z].push_back(timestamp_d);

  return new_z;
}


int LDAweightTOT::sample_z_log(VectorXd &alpha, int &z, int &x,
																		 int &w, int &doc_id)
{
  // removing data is already done in sample_z

  new_z = -1; // debug

	for (int k = 0; k < num_topics; ++k){

		numerator = log( (beta + n_kv(k, w)) *
			(n_dk(doc_id, k) + alpha(k)) );

		denominator = log((double)num_vocab * beta + n_k(k));

		z_prob_vec(k) = numerator - denominator + beta_lg(k);
	}

	sum = logsumexp_Eigen(z_prob_vec, num_topics);
	z_prob_vec = (z_prob_vec.array() - sum).exp();
	new_z = sampler::rcat(z_prob_vec, num_topics);


  // add back data counts
	n_kv(new_z, w) += vocab_weights(w);
	n_k(new_z) += vocab_weights(w);
	n_k_noWeight(new_z) += 1.0;
  n_dk(doc_id, new_z) += 1;

	// Store time stamp
	store_t[new_z].push_back(timestamp_d);

  return new_z;
}



void LDAweightTOT::sample_parameters()
{
	sample_alpha();
	sample_betaparam();
}


void LDAweightTOT::sample_betaparam(){
	// for(int k=0; k<num_topics; k++){
	// 	if(store_t[k].size() == 0)
	// 		continue;
	//
	// 	timestamps_k = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
	// 												(store_t[k].data(), store_t[k].size());		
	//
	// 	beta_mean = timestamps_k.mean();
	// 	beta_var = ( (timestamps_k.array() - beta_mean) * (timestamps_k.array() - beta_mean) ).sum() /
	// 								(store_t[k].size() - 1.0);
	//
	//
	// 	beta_var = 1.0 / (beta_var);  // beta_var reciprocal
	// 	beta_var = ( (beta_mean * (1-beta_mean)) * beta_var - 1.0);
	//	
	// 	beta_params(k, 0) = beta_mean * beta_var;
	// 	beta_params(k, 1) = (1.0 - beta_mean) * beta_var;
	// }

	topic_ids = sampler::shuffled_indexes(num_topics);

	for(int j=0; j<num_topics; j++){
		for(int i=0; i<2; i++){
			k = topic_ids[j];
			current_param = beta_params(k, i);

			start = min_v / (1.0 + min_v);
			end = 1.0;

			previous_p = current_param / (1.0 + current_param);
			slice_ = beta_loglik(k, i) - 2.0 * log(1.0 - previous_p) 
								+ log(unif_rand()); 


			for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++){
				new_p = sampler::slice_uniform(start, end); // <-- using R function above
				beta_params(k, i) = new_p / (1.0 - new_p); // expandp

				newlikelihood = beta_loglik(k, i) - 2.0 * log(1.0 - new_p);

				if (slice_ < newlikelihood){
					store_loglik = newalphallk;
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


double LDAweightTOT::beta_loglik(const int &k, const int &i){
	loglik = 0.0;

  for (int d = 0; d < num_doc; d++){
		loglik += n_dk(d, k) * betapdfln(timestamps(d), beta_params(k,0), beta_params(k,1));
  }

	// Prior
	loglik += gammapdfln(beta_params(k, i), ts_g1, ts_g2);

	return loglik;

}


void LDAweightTOT::sample_alpha()
{

  // start, end, previous_p, new_p, newlikelihood, slice_;
  keep_current_param = alpha;
  topic_ids = sampler::shuffled_indexes(num_topics);
	store_loglik = alpha_loglik();
	newalphallk = 0.0;
	int k;

  for(int i = 0; i < num_topics; i++){
    k = topic_ids[i];
    start = min_v / (1.0 + min_v); // shrinkp
    end = 1.0;
    // end = shrinkp(max_v);
		previous_p = alpha(k) / (1.0 + alpha(k)); // shrinkp
    slice_ = store_loglik - 2.0 * log(1.0 - previous_p) 
						+ log(unif_rand()); // <-- using R random uniform

    for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++){
      new_p = sampler::slice_uniform(start, end); // <-- using R function above
      alpha(k) = new_p / (1.0 - new_p); // expandp

			newalphallk = alpha_loglik();
      newlikelihood = newalphallk - 2.0 * log(1.0 - new_p);

      if (slice_ < newlikelihood){
				store_loglik = newalphallk;
        break;
      } else if (previous_p < new_p){
        end = new_p;
      } else if (new_p < previous_p){
        start = new_p;
      } else {
				Rcpp::stop("Something goes wrong in sample_lambda_slice(). Adjust `A_slice`.");
        alpha(k) = keep_current_param(k);
				break;
      }
    }
  }

	model["alpha"] = alpha;

	// Store alpha
	NumericVector alpha_rvec = alpha_reformat(alpha, num_topics);
	List alpha_iter = model["alpha_iter"];
	alpha_iter.push_back(alpha_rvec);
	model["alpha_iter"] = alpha_iter;
}


double LDAweightTOT::alpha_loglik()
{
  loglik = 0.0;
  fixed_part = 0.0;
  ndk_a = n_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
	alpha_sum_val = alpha.sum();


  fixed_part += lgamma(alpha_sum_val); // first term numerator
  for(int k = 0; k < num_topics; k++){
    fixed_part -= lgamma(alpha(k)); // first term denominator
    // Add prior
		loglik += gammapdfln(alpha(k), eta_1_regular, eta_2_regular);

  }

  for(int d = 0; d < num_doc; d++){
    loglik += fixed_part;
    // second term numerator
    for(int k = 0; k < num_topics; k++){
      loglik += mylgamma(ndk_a(d,k));
    }
    // second term denominator
    loglik -= mylgamma(doc_each_len[d] + alpha_sum_val);

  }
  return loglik;
}


double LDAweightTOT::loglik_total()
{
  double loglik = 0.0;
  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += mylgamma(beta + n_kv(k, v) / vocab_weights(v) ) - mylgamma(beta);
    }
    // word normalization
    loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_k_noWeight(k) );

		// Rcout << (double)n_x0_k(k) << " / " << (double)n_x1_k(k) << std::endl; // debug
  }
  // z and time stamps
	fixed_part = alpha.sum();
  for (int d = 0; d < num_doc; d++){
    loglik += mylgamma( fixed_part ) - mylgamma( doc_each_len[d] + fixed_part );
    for (int k = 0; k < num_topics; k++){
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );

			// time stamps
			loglik += n_dk(d, k) * betapdfln(timestamps(d), beta_params(k,0), beta_params(k,1));
    }
  }

	// Prior
	for(int k=0; k<num_topics; k++){
		for(int i=0; i<2; i++){
			loglik += gammapdfln(beta_params(k, i), ts_g1, ts_g2);
		}	
	}

  return loglik;
}


void LDAweightTOT::verbose_special(int &r_index){
	// Store beta param

	Rcpp::NumericMatrix tot_beta_R = Rcpp::wrap(beta_params);
	List tot_beta = model["tot_beta"];
	tot_beta.push_back(tot_beta_R);
	model["tot_beta"] = tot_beta;
}


