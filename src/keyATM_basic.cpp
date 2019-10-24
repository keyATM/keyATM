#include "keyATM_basic.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


void keyATMbasic::read_data_specific()
{
  nv_alpha = priors_list["alpha"];
  alpha = Rcpp::as<Eigen::VectorXd>(nv_alpha);

  prior_gamma = MatrixXd::Zero(num_topics, 2);
  NumericMatrix RMatrix = priors_list["gamma"];
  prior_gamma = Rcpp::as<Eigen::MatrixXd>(RMatrix);
  beta_s = priors_list["beta_s"];

  estimate_alpha = options_list["estimate_alpha"];
  if (estimate_alpha == 0) {
    store_alpha = 0;
  } else {
    store_alpha = 1;
  }
}


void keyATMbasic::initialize_specific()
{
  // No additional initialization
}

void keyATMbasic::iteration_single(int &it)
{ // Single iteration

  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  for (int ii = 0; ii < num_doc; ii++){
    doc_id_ = doc_indexes[ii];
    doc_x = X[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];
    
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
  sample_parameters(it);

}

void keyATMbasic::sample_parameters(int &it)
{
  if (estimate_alpha)
    sample_alpha();

  // Store alpha
  if (store_alpha){
    int r_index = it + 1;
    if (r_index % thinning == 0 || r_index == 1 || r_index == iter) {
      NumericVector alpha_rvec = alpha_reformat(alpha, num_topics);
      List alpha_iter = stored_values["alpha_iter"];
      alpha_iter.push_back(alpha_rvec);
      stored_values["alpha_iter"] = alpha_iter;  
    }
  }
}


void keyATMbasic::sample_alpha()
{

  // start, end, previous_p, new_p, newlikelihood, slice_;
  keep_current_param = alpha;
  topic_ids = sampler::shuffled_indexes(num_topics);
  store_loglik = alpha_loglik();
  newalphallk = 0.0;
  int k;

  for (int i = 0; i < num_topics; i++) {
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
}


double keyATMbasic::alpha_loglik()
{
  loglik = 0.0;
  
  fixed_part = 0.0;
  ndk_a = n_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
  alpha_sum_val = alpha.sum();
  
  
  fixed_part += mylgamma(alpha_sum_val); // first term numerator
  for (int k = 0; k < num_topics; k++) {
    fixed_part -= mylgamma(alpha(k)); // first term denominator
    // Add prior
    if (k < keyword_k) {
      loglik += gammapdfln(alpha(k), eta_1, eta_2);
    } else {
      loglik += gammapdfln(alpha(k), eta_1_regular, eta_2_regular);
    }
  
  }
  for (int d = 0; d < num_doc; d++) {
    loglik += fixed_part;
    // second term numerator
    for (int k = 0; k < num_topics; k++) {
      loglik += mylgamma(ndk_a(d,k));
    }
    // second term denominator
    loglik -= mylgamma(doc_each_len[d] + alpha_sum_val);
  
  }

  return loglik;
}


double keyATMbasic::loglik_total()
{
  double loglik = 0.0;


  for (int k = 0; k < num_topics; k++){
    for (int v = 0; v < num_vocab; v++){ // word
      loglik += mylgamma(beta + n_x0_kv(k, v) / vocab_weights(v) ) - mylgamma(beta);
      // loglik += mylgamma(beta_s + n_x1_kv.coeffRef(k, v) / vocab_weights(v) ) - mylgamma(beta_s);
    }



    // word normalization
    loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_x0_k_noWeight(k) );

    if (k < keyword_k) {
      // For keyword topics

      // n_x1_kv
      for (SparseMatrix<double,RowMajor>::InnerIterator it(n_x1_kv, k); it; ++it){
        loglik += mylgamma(beta_s + it.value() / vocab_weights(it.index()) ) - mylgamma(beta_s);
      }
      loglik += mylgamma( beta_s * (double)keywords_num[k] ) - mylgamma(beta_s * (double)keywords_num[k] + n_x1_k_noWeight(k) );

      // Normalization
      loglik += mylgamma( prior_gamma(k, 0) + prior_gamma(k, 1)) - mylgamma( prior_gamma(k, 0)) - mylgamma( prior_gamma(k, 1));

      // x
      loglik += mylgamma( n_x0_k_noWeight(k) + prior_gamma(k, 1) ) 
                -  mylgamma(n_x1_k_noWeight(k) + prior_gamma(k, 0) + n_x0_k_noWeight(k) + prior_gamma(k, 1))
                + mylgamma( n_x1_k_noWeight(k) + prior_gamma(k, 0) );  
    }
  }
  // z
  fixed_part = alpha.sum();
  for (int d = 0; d < num_doc; d++){
    loglik += mylgamma( fixed_part ) - mylgamma( doc_each_len[d] + fixed_part );
    for (int k = 0; k < num_topics; k++){
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }
  }

  return loglik;
}

