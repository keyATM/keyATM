#include "keyATM_covPG.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */


void keyATMcovPG::read_data_specific()
{
  // Covariate
  model_settings = model["model_settings"];
  NumericMatrix C_r = model_settings["covariates_data_use"];
  num_cov = C_r.cols();

  // Access to PG related parameters
  PG_params = model_settings["PG_params"];
}

void keyATMcovPG::initialize_specific()
{
  theta = MatrixXd::Zero(num_doc, num_topics);
}


void keyATMcovPG::iteration_single(int it)
{ // Single iteration
  int doc_id_;
  int doc_length;
  int w_, z_, s_;
  int new_z, new_s;
  int w_position;


  sample_parameters(it);

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
    
      new_z = sample_z_PG(z_, s_, w_, doc_id_);
      doc_z[w_position] = new_z;
    
      if (keywords[new_z].find(w_) == keywords[new_z].end())	
        continue;
  
      z_ = doc_z[w_position]; // use updated z
      new_s = sample_s(z_, s_, w_, doc_id_);
      doc_s[w_position] = new_s;
    }
    
    Z[doc_id_] = doc_z;
    S[doc_id_] = doc_s;
  }
}


void keyATMcovPG::sample_parameters(int it)
{
  sample_PG(it);

  // Store theta 
  int r_index = it + 1;
  if (store_theta) {
    if (r_index % thinning == 0 || r_index == 1 || r_index == iter) {
      Rcpp::NumericMatrix theta_R = Rcpp::wrap(theta);
      List theta_iter = stored_values["theta_PG"];
      theta_iter.push_back(theta_R);
      stored_values["theta_PG"] = theta_iter;
    }
  }

  if (r_index == iter) {
    PG_params["theta_last"] = Rcpp::wrap(theta);
    model_settings["PG_params"] = PG_params;
  }
}


void keyATMcovPG::sample_PG(int it)
{
  // multiPGreg function
  Environment pkg = Environment::namespace_env("keyATM");
  Function PGreg_Rfun = pkg["multiPGreg"];
  NumericMatrix C_r = model_settings["covariates_data_use"];
  NumericMatrix Y_r = Rcpp::wrap(n_dk);

  int r_index = it + 1;
  int store_lambda = 0;
  if (r_index % thinning == 0 || r_index == 1 || r_index == iter) {
    store_lambda = 1;
  }

  PG_params = PGreg_Rfun(Y_r, C_r, num_topics, PG_params, 1, store_lambda);
  NumericMatrix theta_tilda_r = PG_params["theta_tilda"];
  utils::calc_PGtheta(theta_tilda_r, theta, num_doc, num_topics);  // update theta
}


int keyATMcovPG::sample_z_PG(int z, int s, int w, int doc_id)
{
  int new_z;
  double numerator, denominator;
  double sum;

  // remove data
  if (s == 0) {
    n_s0_kv(z, w) -= vocab_weights(w);
    n_s0_k(z) -= vocab_weights(w);
  } else if (s == 1) {
    n_s1_kv.coeffRef(z, w) -= vocab_weights(w);
    n_s1_k(z) -= vocab_weights(w);
  } else {
    Rcerr << "Error at sample_z, remove" << std::endl;
  }

  n_dk(doc_id, z) -= vocab_weights(w);
  n_dk_noWeight(doc_id, z) -= 1.0;

  new_z = -1; // debug
  if (s == 0) {
    for (int k = 0; k < num_topics; ++k) {

      numerator = (beta + n_s0_kv(k, w)) *
        (n_s0_k(k) + prior_gamma(k, 1)) *
        theta(doc_id, k);

      denominator = (Vbeta + n_s0_k(k)) *
        (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1));

      z_prob_vec(k) = numerator / denominator;
    }

    sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  } else {
    for (int k = 0; k < num_topics; ++k) {
      if (keywords[k].find(w) == keywords[k].end()) {
        z_prob_vec(k) = 0.0;
        continue;
      } else { 
        numerator = (beta_s + n_s1_kv.coeffRef(k, w)) *
          (n_s1_k(k) + prior_gamma(k, 0)) *
          theta(doc_id, k);
        denominator = (Lbeta_sk(k) + n_s1_k(k) ) *
          (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1));

        z_prob_vec(k) = numerator / denominator;
      }
    }

    sum = z_prob_vec.sum();
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  }

  // add back data counts
  if (s == 0) {
    n_s0_kv(new_z, w) += vocab_weights(w);
    n_s0_k(new_z) += vocab_weights(w);
  } else if (s == 1) {
    n_s1_kv.coeffRef(new_z, w) += vocab_weights(w);
    n_s1_k(new_z) += vocab_weights(w);
  } else {
    Rcerr << "Error at sample_z, add" << std::endl;
  }
  n_dk(doc_id, new_z) += vocab_weights(w);
  n_dk_noWeight(doc_id, new_z) += 1.0;

  return new_z;
}


double keyATMcovPG::loglik_total()
{
  double loglik = 0.0;
  for (int k = 0; k < num_topics; ++k) {
    for (int v = 0; v < num_vocab; ++v) { // word
      loglik += mylgamma(beta + n_s0_kv(k, v)) - mylgamma(beta);
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
                + mylgamma( n_s1_k(k) + prior_gamma(k, 0) );  
    }
  }


  // z
  for (int d = 0; d < num_doc; ++d) {
    for (int k = 0; k < num_topics; ++k) {
      loglik += log(theta(d, k)) * n_dk(d, k);
    }
  }

  return loglik;
}



