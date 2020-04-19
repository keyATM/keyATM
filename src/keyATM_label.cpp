#include "keyATM_label.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


void keyATMlabel::read_data_specific()
{
  nv_alpha = priors_list["alpha"];
  alpha = Rcpp::as<Eigen::VectorXd>(nv_alpha);

  // read label data
  model_settings = model["model_settings"];
  label_vec = model_settings["labels"];

  estimate_alpha = options_list["estimate_alpha"];
  if (estimate_alpha == 0) {
    store_alpha = 0;
  } else {
    store_alpha = 1;
  }
}


void keyATMlabel::initialize_specific()
{
  // Initialize label information
  label_dk = MatrixXd::Zero(num_doc, num_topics);

  // when use non-weighted length
  for (int i = 0; i < num_doc; i++) {
    doc_label = label_vec[i];
  // if the label is less than zero, it means label is missing
    if (doc_label >= 0) {
        label_dk(i, doc_label) = doc_each_len_weighted[i]; // use non-log weighted doc-length
    }
  }
  // Alpha to store during the iteration
  Alpha = MatrixXd::Zero(num_doc, num_topics);
}

void keyATMlabel::iteration_single(int it)
{ // Single iteration
  int doc_id_;
  int doc_length;
  int w_, z_, s_;
  int new_z, new_s;
  int w_position;

  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle
  Alpha = label_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting

  for (int ii = 0; ii < num_doc; ii++) {
    doc_id_ = doc_indexes[ii];
    doc_s = S[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_];
    doc_length = doc_each_len[doc_id_];

    alpha = Alpha.row(doc_id_).transpose(); // chooose document specific alpha
    token_indexes = sampler::shuffled_indexes(doc_length); //shuffle
    
    // Iterate each word in the document
    for (int jj = 0; jj < doc_length; jj++) {
      w_position = token_indexes[jj];
      s_ = doc_s[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position];
    
      new_z = (use_labels) ? sample_z_label(alpha, z_, s_, w_, doc_id_) : sample_z(alpha, z_, s_, w_, doc_id_);
      doc_z[w_position] = new_z;
    
      if (keywords[new_z].find(w_) == keywords[new_z].end())	
        continue;
  
      z_ = doc_z[w_position]; // use updated z
      new_s = (use_labels) ? sample_s_label(alpha, z_, s_, w_, doc_id_) : sample_s(alpha, z_, s_, w_, doc_id_);
      doc_s[w_position] = new_s;
    }
    
    Z[doc_id_] = doc_z;
    S[doc_id_] = doc_s;
  }
  sample_parameters(it);
}


void keyATMlabel::sample_parameters(int it)
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


void keyATMlabel::sample_alpha()
{

  keep_current_param = alpha;
  topic_ids = sampler::shuffled_indexes(num_topics);

  newalphallk = 0.0;
  int k;

  for (int i = 0; i < num_topics; i++) {
    k = topic_ids[i];
    store_loglik = alpha_loglik_label(k);
    start = min_v / (1.0 + min_v); // shrinkp
    end = 1.0;
    // end = shrinkp(max_v);
    previous_p = alpha(k) / (1.0 + alpha(k)); // shrinkp
    slice_ = store_loglik - 2.0 * log(1.0 - previous_p) 
            + log(unif_rand()); // <-- using R random uniform

    for (int shrink_time = 0; shrink_time < max_shrink_time; shrink_time++) {
      new_p = sampler::slice_uniform(start, end); // <-- using R function above
      alpha(k) = new_p / (1.0 - new_p); // expandp

      newalphallk = alpha_loglik_label(k);
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


double keyATMlabel::alpha_loglik_label(int k)
{
  loglik = 0.0;
  
  // fixed_part = 0.0;
  Alpha = label_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
  ndk_a = n_dk + Alpha; // Adding two matrices

  Alpha_sum_vec = Alpha.rowwise().sum(); 

  // No fixed part this time 
  // Add prior
  if (k < keyword_k) {
    loglik += gammapdfln(alpha(k), eta_1, eta_2);
  } else {
    loglik += gammapdfln(alpha(k), eta_1_regular, eta_2_regular);
  }

  for (int d = 0; d < num_doc; d++) {
    // loglik += fixed_part;
    loglik += mylgamma(Alpha_sum_vec(d)); // first term numerator
    loglik -= mylgamma(Alpha(d, k)); // first term denominator

    // second term numerator
    loglik += mylgamma(ndk_a(d,k));

    // second term denominator
    loglik -= mylgamma(doc_each_len_weighted[d] + Alpha_sum_vec(d));
  }

  return loglik;
}


double keyATMlabel::loglik_total()
{
  double loglik = 0.0;

  for (int k = 0; k < num_topics; k++) {
    for (int v = 0; v < num_vocab; v++) { // word
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
  Alpha = label_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
  ndk_a = n_dk + Alpha; // Adding two matrices
  Alpha_sum_vec = Alpha.rowwise().sum();

  for (int d = 0; d < num_doc; d++) {
    loglik += mylgamma( Alpha_sum_vec(d) ) - mylgamma( doc_each_len_weighted[d] + Alpha_sum_vec(d) );

    for (int k = 0; k < num_topics; k++) {
      loglik += mylgamma( ndk_a(d,k) ) - mylgamma( Alpha(d, k) );
    }
  }

  return loglik;
}


double keyATMlabel::loglik_total_label()
{
  double loglik = 0.0;

  for (int k = 0; k < num_topics; k++) {
    for (int v = 0; v < num_vocab; v++) { // word
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
  Alpha = label_dk.rowwise() + alpha.transpose(); // Use Eigen Broadcasting
  ndk_a = n_dk + Alpha; // Adding two matrices
  Alpha_sum_vec = Alpha.rowwise().sum();

  for (int d = 0; d < num_doc; d++) {
    loglik += mylgamma( Alpha_sum_vec(d) ) - mylgamma( doc_each_len_weighted[d] + Alpha_sum_vec(d) );

    for (int k = 0; k < num_topics; k++) {
      loglik += mylgamma( ndk_a(d,k) ) - mylgamma( Alpha(d, k) );
    }
  }

  return loglik;
}
