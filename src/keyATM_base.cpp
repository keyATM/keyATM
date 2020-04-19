#include "keyATM_base.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


void keyATMbase::read_data_specific()
{
  nv_alpha = priors_list["alpha"];
  alpha = Rcpp::as<Eigen::VectorXd>(nv_alpha);

  estimate_alpha = options_list["estimate_alpha"];
  if (estimate_alpha == 0) {
    store_alpha = 0;
  } else {
    store_alpha = 1;
  }
}


void keyATMbase::initialize_specific()
{
  // No additional initialization
  // This part is used when there is a model specific need
  // to initialize variables.
}

void keyATMbase::iteration_single(int it)
{ // Single iteration
  int doc_id_;
  int doc_length;
  int w_, z_, s_;
  int new_z, new_s;
  int w_position;

  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  for (int ii = 0; ii < num_doc; ++ii) {
    doc_id_ = doc_indexes[ii];
    doc_s = S[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];
    
    token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
    
    // Iterate each word in the document
    for (int jj = 0; jj < doc_length; ++jj) {
      w_position = token_indexes[jj];
      s_ = doc_s[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position];
    
      new_z = sample_z(alpha, z_, s_, w_, doc_id_);
      doc_z[w_position] = new_z;
    
      if (keywords[new_z].find(w_) == keywords[new_z].end())	
        continue;
  
      z_ = doc_z[w_position]; // use updated z
      new_s = sample_s(alpha, z_, s_, w_, doc_id_);
      doc_s[w_position] = new_s;
    }
    
    Z[doc_id_] = doc_z;
    S[doc_id_] = doc_s;
  }
  sample_parameters(it);

}

void keyATMbase::sample_parameters(int it)
{
  if (estimate_alpha)
    sample_alpha();

  // Store alpha
  if (store_alpha) {
    int r_index = it + 1;
    if (r_index % thinning == 0 || r_index == 1 || r_index == iter) {
      NumericVector alpha_rvec = alpha_reformat(alpha, num_topics);
      List alpha_iter = stored_values["alpha_iter"];
      alpha_iter.push_back(alpha_rvec);
      stored_values["alpha_iter"] = alpha_iter;  
    }
  }
}


void keyATMbase::sample_alpha()
{

  double start, end, previous_p, new_p, newlikelihood, slice_;
  keep_current_param = alpha;
  topic_ids = sampler::shuffled_indexes(num_topics);
  newalphallk = 0.0;
  int k;

  for (int i = 0; i < num_topics; ++i) {
    k = topic_ids[i];
    store_loglik = alpha_loglik(k);
    start = min_v ; // shrinked with shrinkp()
    end = max_v;  // shrinked with shrinkp()

    previous_p = alpha(k) / (1.0 + alpha(k)); // shrinkp
    slice_ = store_loglik - 2.0 * log(1.0 - previous_p) 
            + log(unif_rand()); // <-- using R random uniform

    for (int shrink_time = 0; shrink_time < max_shrink_time; ++shrink_time) {
      new_p = sampler::slice_uniform(start, end); // <-- using R function above
      alpha(k) = new_p / (1.0 - new_p); // expandp

      newalphallk = alpha_loglik(k);
      newlikelihood = newalphallk - 2.0 * log(1.0 - new_p);

      if (slice_ < newlikelihood) {
        break;
      } else if (previous_p < new_p) {
        end = new_p;
      } else if (new_p < previous_p) {
        start = new_p;
      } else {
        Rcpp::stop("Something goes wrong in sample_lambda_slice().");
        alpha(k) = keep_current_param(k);
        break;
      }
    }
  }
}


double keyATMbase::alpha_loglik(int k)
{
  double loglik = 0.0;
  double fixed_part = 0.0;

  ndk_a = n_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
  double alpha_sum_val = alpha.sum();
  
  
  fixed_part += mylgamma(alpha_sum_val); // first term numerator
  fixed_part -= mylgamma(alpha(k)); // first term denominator
  // Add prior
  if (k < keyword_k) {
    loglik += gammapdfln(alpha(k), eta_1, eta_2);
  } else {
    loglik += gammapdfln(alpha(k), eta_1_regular, eta_2_regular);
  }

  for (int d = 0; d < num_doc; ++d) {
    loglik += fixed_part;

    // second term numerator
    loglik += mylgamma(ndk_a(d,k));

    // second term denominator
    loglik -= mylgamma(doc_each_len_weighted[d] + alpha_sum_val);
  }

  return loglik;
}


double keyATMbase::loglik_total()
{
  double loglik = 0.0;
  double fixed_part = 0.0;

  for (int k = 0; k < num_topics; ++k) {
    for (int v = 0; v < num_vocab; ++v) { // word
      loglik += mylgamma(beta + n_s0_kv(k, v) ) - mylgamma(beta);
    }

    // word normalization
    loglik += mylgamma( beta * (double)num_vocab ) - mylgamma(beta * (double)num_vocab + n_s0_k(k) );

    if (k < keyword_k) {
      // For keyword topics

      // n_s1_kv
      for (SparseMatrix<double,RowMajor>::InnerIterator it(n_s1_kv, k); it; ++it) {
        loglik += mylgamma(beta_s + it.value()) - mylgamma(beta_s);
      }
      loglik += mylgamma( beta_s * (double)keywords_num[k] ) - mylgamma(beta_s * (double)keywords_num[k] + n_s1_k(k) );

      // Normalization
      loglik += mylgamma( prior_gamma(k, 0) + prior_gamma(k, 1)) - mylgamma( prior_gamma(k, 0)) - mylgamma( prior_gamma(k, 1));

      // s
      loglik += mylgamma( n_s0_k(k) + prior_gamma(k, 1) ) 
                - mylgamma(n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1))
                + mylgamma(n_s1_k(k) + prior_gamma(k, 0));  
    }
  }

  // z
  fixed_part = alpha.sum();
  for (int d = 0; d < num_doc; ++d) {
    loglik += mylgamma( fixed_part ) - mylgamma( doc_each_len_weighted[d] + fixed_part );

    for (int k = 0; k < num_topics; ++k) {
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }
  }

  return loglik;
}

double keyATMbase::loglik_total_label()
{
  double loglik = 0.0;
  double fixed_part = 0.0;

  for (int k = 0; k < num_topics; ++k) {
    for (int v = 0; v < num_vocab; ++v) { // word
      loglik += mylgamma(beta_s0kv(k, v) + n_s0_kv(k, v) ) - mylgamma(beta_s0kv(k, v));
    }

    // word normalization
    loglik += mylgamma( Vbeta_k(k) ) - mylgamma(Vbeta_k(k) + n_s0_k(k) );

    if (k < keyword_k) {
      // For keyword topics

      // n_s1_kv
      for (SparseMatrix<double,RowMajor>::InnerIterator it(n_s1_kv, k); it; ++it) {
        loglik += mylgamma(beta_s1kv.coeffRef(k, it.index()) + it.value()) - mylgamma(beta_s1kv.coeffRef(k, it.index()));
      }
      loglik += mylgamma( Lbeta_sk(k) ) - mylgamma(Lbeta_sk(k) + n_s1_k(k) );

      // Normalization
      loglik += mylgamma( prior_gamma(k, 0) + prior_gamma(k, 1)) - mylgamma( prior_gamma(k, 0)) - mylgamma( prior_gamma(k, 1));

      // s
      loglik += mylgamma( n_s0_k(k) + prior_gamma(k, 1) ) 
                - mylgamma(n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1))
                + mylgamma(n_s1_k(k) + prior_gamma(k, 0));  
    }
  }

  // z
  fixed_part = alpha.sum();
  for (int d = 0; d < num_doc; ++d) {
    loglik += mylgamma( fixed_part ) - mylgamma( doc_each_len_weighted[d] + fixed_part );

    for (int k = 0; k < num_topics; ++k) {
      loglik += mylgamma( n_dk(d,k) + alpha(k) ) - mylgamma( alpha(k) );
    }
  }

  return loglik;
}
