#include "keyATM_tot.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


keyATMtot::keyATMtot(List model_, const int iter_, const int output_per_) :
	keyATMbase(model_, iter_, output_per_) // pass to parent!
{

	// Constructor

	read_data();
	initialize();
	iteration();
}


void keyATMtot::read_data_specific()
{
	alpha = Rcpp::as<Eigen::VectorXd>(nv_alpha);

	// Read time stamp
	timestamps = Rcpp::as<Eigen::VectorXd>(model["timestamps"]);
}


void keyATMtot::initialize_specific()
{
	// Store parameters for time Beta
	beta_params = MatrixXd::Constant(num_topics, 2, 0.5);
	// beta_tg_num = VectorXd::Zero(num_topics);
	// beta_tg_denom = VectorXd::Zero(num_topics);
	// beta_lg_num = VectorXd::Zero(num_topics);
	// beta_lg_denom = VectorXd::Zero(num_topics);

	beta_tg = VectorXd::Zero(num_topics);
	beta_lg = VectorXd::Zero(num_topics);
	beta_tg_base = VectorXd::Zero(num_topics);
	beta_lg_base = VectorXd::Zero(num_topics);

	// Store time stamps
	for(int k=0; k<num_topics; k++){
		vector<double> temp;
		store_t.push_back(temp);
	}

	Rcerr << "Initialization specific ends." << endl;  // debug

}

void keyATMtot::iteration_single(int &it)
{ // Single iteration

	doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle
	Rcerr << "Iteration starts" << endl;  // debug

	// Clear temporary time stamp
	for(int k=0; k<num_topics; k++){
		store_t[k].clear();
	}

	// Apply gamma
	for(int k=0; k<num_topics; k++){ 
		beta_a = beta_params(k, 0);
		beta_b = beta_params(k, 1);	

		// Rcerr << "Beta a/b: " << beta_a << "  /  " << beta_b << endl;  // debug
	
		// Log version
		beta_lg(k) = beta_lg_base(k) + mylgamma(beta_a + beta_b) - (mylgamma(beta_a) + mylgamma(beta_b));	

		// Normal version
		beta_tg(k) = beta_tg_base(k) * tgamma(beta_a + beta_b) / (tgamma(beta_a) * tgamma(beta_b));


		// // Log version
		// beta_lg_num(k) = mylgamma(beta_a + beta_b);
		// beta_lg_denom(k) = mylgamma(beta_a) + mylgamma(beta_b);
		//
		// // Normal version
		// beta_tg_num(k) = tgamma(beta_a + beta_b);
		// beta_tg_denom(k) = tgamma(beta_a) * tgamma(beta_b);
		//
		//
		// Rcerr <<  beta_lg_num(k) << " " << beta_lg_denom(k)   << " " <<  beta_tg_num(k)  << " " << beta_tg_denom(k) << endl;
		//
		//
		// 	if(beta_lg_num(k) < numeric_limits<double>::min() | 
		// 			beta_lg_num(k) > numeric_limits<double>::max()){
		// 		Rcerr << "Overflow/underflow beta_lg_num " << beta_lg_num(k)  << endl;  // debug
		// 	}
		//
		// 	if(beta_lg_denom(k) < numeric_limits<double>::min() | 
		// 			beta_lg_denom(k) > numeric_limits<double>::max()){
		// 		Rcerr << "Overflow/underflow beta_lg_denom " << beta_lg_denom(k) << endl;  // debug
		// 	}
		//
		// 	if(beta_tg_num(k) < numeric_limits<double>::min() | 
		// 			beta_tg_num(k) > numeric_limits<double>::max()){
		// 		Rcerr << "Overflow/underflow beta_tg_num " << beta_tg_num(k) << endl;  // debug
		// 	}
		//
		// 	if(beta_tg_denom(k) < numeric_limits<double>::min() | 
		// 			beta_tg_denom(k) > numeric_limits<double>::max()){
		// 		Rcerr << "Overflow/underflow beta_tg_denom " << beta_tg_denom(k) << endl;  // debug
		// 	}
	}

	

	// Sampling
	for (int ii = 0; ii < num_doc; ii++){
		doc_id_ = doc_indexes[ii];
		doc_x = X[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
		doc_length = doc_each_len[doc_id_];
		
		token_indexes = sampler::shuffled_indexes(doc_length); //shuffle

		// Prepare beta_a and beta_b for sampling
		timestamp_d = timestamps(doc_id_);

		for(int k=0; k<num_topics; k++){ 
			beta_a = beta_params(k, 0);
			beta_b = beta_params(k, 1);
		
			// for log version
			beta_lg(k) += (beta_a - 1.0) * log(1.0 - timestamp_d) + (beta_b - 1.0) * log(timestamp_d);
		
			// For normal version
			beta_tg(k) *=  pow(1.0 - timestamp_d, beta_a - 1.0) * pow(timestamp_d, beta_b - 1.0);
			
		}
		
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

	cout << beta_params << endl;

	// if(floor(iter/2) < it)

}


int keyATMtot::sample_z_log(VectorXd &alpha, int &z, int &x,
																		 int &w, int &doc_id)
{
  // removing data is already done in sample_z
 //  if (x == 0){
 //    n_x0_kv(z, w) -= vocab_weights(w);
 //    n_x0_k(z) -= vocab_weights(w);
 //    n_x0_k_noWeight(z) -= 1.0;
 //  } else if (x==1) {
 //    n_x1_kv.coeffRef(z, w) -= vocab_weights(w);
 //    n_x1_k(z) -= vocab_weights(w);
 //    n_x1_k_noWeight(z) -= 1.0;
 //  } else {
	// 	Rcerr << "Error at sample_z, remove" << std::endl;
	// }
	//
 //  n_dk(doc_id, z) -= 1;

  new_z = -1; // debug
  if (x == 0){
    for (int k = 0; k < num_topics; ++k){

      numerator = log( (beta + n_x0_kv(k, w)) *
        (n_x0_k(k) + gamma_2) *
        (n_dk(doc_id, k) + alpha(k)) );
				// beta_lg_num(k);

      denominator = log( ((double)num_vocab * beta + n_x0_k(k)) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2) );
				// beta_lg_denom(k);

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
						// beta_lg_num(k);
      }
      denominator = log( ((double)seed_num[k] * beta_s + n_x1_k(k) ) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2) );
				// beta_lg_denom(k);

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


int keyATMtot::sample_z(VectorXd &alpha, int &z, int &x,
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

	//debug
	// new_z = z;
	// goto UPDATE;
	

  if (x == 0){
    for (int k = 0; k < num_topics; ++k){
      numerator = (beta + n_x0_kv(k, w)) *
        (n_x0_k(k) + gamma_2) *
        (n_dk(doc_id, k) + alpha(k));
				// beta_tg_num(k);

      denominator = ((double)num_vocab * beta + n_x0_k(k)) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2);
				// beta_tg_denom(k);


			check_frac = numerator / denominator * beta_tg(k);

			if(check_frac < numeric_limits<double>::min() | 
					check_frac > numeric_limits<double>::max()){
				// Rcerr << "Overflow/underflow" << endl;  // debug
				// new_z = z;
				// goto UPDATE;
				return sample_z_log(alpha, z, x, w, doc_id);
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
						// beta_tg_num(k);
      }
      denominator = ((double)seed_num[k] * beta_s + n_x1_k(k) ) *
        (n_x1_k(k) + gamma_1 + n_x0_k(k) + gamma_2);
				// beta_tg_denom(k);


			check_frac = numerator / denominator * beta_tg(k);

			if(check_frac < numeric_limits<double>::min() | 
					check_frac > numeric_limits<double>::max()){
				// Rcerr << "Overflow/underflow" << endl;  // debug
				// new_z = z;
				// goto UPDATE;
				return sample_z_log(alpha, z, x, w, doc_id);
			}

			z_prob_vec(k) = check_frac;
    }


		sum = z_prob_vec.sum();
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  }

  // add back data counts

	UPDATE:
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



void keyATMtot::sample_parameters()
{
	sample_alpha();
	sample_betaparam();
}


void keyATMtot::sample_betaparam(){
	Rcerr << "Sample beta_params" << endl;

	for(int k=0; k<num_topics; k++){
		timestamps_k = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>
													(store_t[k].data(), store_t[k].size());		
	
		beta_mean = timestamps_k.mean();
		beta_var = ( (timestamps_k.array() - beta_mean) * (timestamps_k.array() - beta_mean) ).sum() /
									(store_t[k].size() - 1.0);
	
	
		beta_var = 1.0 / (beta_var);  // beta_var reciprocal
		beta_var = ( (beta_mean * (1-beta_mean)) * beta_var - 1.0);
		
		// cout << beta_mean * beta_var << "/" << (1.0 - beta_mean) * beta_var << endl;  // debug
		beta_params(k, 0) = beta_mean * beta_var;
		beta_params(k, 1) = (1.0 - beta_mean) * beta_var;
	}
}


void keyATMtot::sample_alpha()
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


void keyATMtot::verbose_special(int &r_index){
	// Store beta param

	Rcpp::NumericMatrix tot_beta_R = Rcpp::wrap(beta_params);
	List tot_beta = model["tot_beta"];
	tot_beta.push_back(tot_beta_R);
	model["tot_beta"] = tot_beta;
}


double keyATMtot::alpha_loglik()
{
  loglik = 0.0;
	
  fixed_part = 0.0;
  ndk_a = n_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
	alpha_sum_val = alpha.sum();
	
	
  fixed_part += mylgamma(alpha_sum_val); // first term numerator
  for(int k = 0; k < num_topics; k++){
    fixed_part -= mylgamma(alpha(k)); // first term denominator
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
      loglik += mylgamma(ndk_a(d,k));
    }
    // second term denominator
    loglik -= mylgamma(doc_each_len[d] + alpha_sum_val);
	
  }

  return loglik;
}


double keyATMtot::loglik_total()
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

		// Rcout << (double)n_x0_k(k) << " / " << (double)n_x1_k(k) << std::endl; // debug

    // x normalization
    loglik += mylgamma(gamma_1 + gamma_2) - mylgamma(gamma_1) - mylgamma(gamma_2);

  }
  // z
  for (int d = 0; d < num_doc; d++){
    loglik += mylgamma( alpha.sum() ) - mylgamma( doc_each_len[d] + alpha.sum() );
    for (int k = 0; k < num_topics; k++){
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );


			// time stamps
			loglik += n_dk(d, k) * betapdfln(timestamps(d), beta_params(k,0), beta_params(k,1));
    }
  }
  return loglik;
}


