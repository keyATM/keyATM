#include "keyATM_multicov.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

# define PI_V   3.14159265358979323846  /* pi */

void keyATMmulticov::read_data_specific()
{
  num_corpus = model_settings["num_corpus"]; //Nr. of corpora
  corpus_id = model_settings["corpus_id"]; // Vector of ints of corpus id for each document
  num_doc_all = model_settings["docs_per_corpus"]; //Vector of ints, nr. of documents per corpus
  doc_lookup = stored_values["doc_lookup"];
  doc_lookup_global = stored_values["doc_lookup_global"];

  covariate_data = model_settings["covariates_data_use"];

  C = model["C"];

  X.resize(num_corpus);
  num_cov_all = VectorXd::Zero(num_corpus);

  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    NumericMatrix covariate_data_corpus = covariate_data[std::to_string(corpus + 1)];
    X[corpus] = Rcpp::as<Eigen::MatrixXd>(covariate_data_corpus);
    num_cov_all[corpus] = X[corpus].cols();
  }

  NumericMatrix RMatrix = priors_list["omega"];
  prior_omega = Rcpp::as<Eigen::MatrixXd>(RMatrix);

  beta_c = priors_list["beta_c"];
  // Slice Sampling
  val_min = model_settings["slice_min"];
  val_min = shrink(val_min, slice_A);

  val_max = model_settings["slice_max"];
  val_max = shrink(val_max, slice_A);

  mu_all = VectorXd::Zero(num_corpus);
  sigma_all = VectorXd::Zero(num_corpus);
  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    sigma_all[corpus] = 1.0;
  }
}

void keyATMmulticov::create_sufficient_stats()
{
  Vbeta_c = double(num_vocab) * beta_c;
  Vbeta_s = double(num_vocab) * beta_s;

  int s_, z_, w_, c_, corpus_, doc_len_;
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  n_dk_all.resize(num_corpus);

  n_s0_c0_kv_all.resize(num_corpus);
  
  n_s0_c0_k_all = MatrixXd::Zero(num_topics, num_corpus);
  n_s0_c1_k = VectorXd::Zero(num_topics);
  n_s0_c1_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_c0_k = VectorXd::Zero(num_topics);
  n_c1_k = VectorXd::Zero(num_topics);
  n_c0_k_all = MatrixXd::Zero(num_topics, num_corpus);
  n_c1_k_all = MatrixXd::Zero(num_topics, num_corpus);
  
  n_s0_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_s1_kv_multi.resize(num_topics, num_vocab);
  n_s0_k = VectorXd::Zero(num_topics);
  n_s1_k = VectorXd::Zero(num_topics);

  vocab_weights_corpus = MatrixXd::Constant(num_corpus, num_vocab, 1.0);
  total_words_corpus.resize(num_corpus);

  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    n_s0_c0_kv_all[corpus] = MatrixXd::Zero(num_topics, num_vocab);
    n_dk_all[corpus] = MatrixXd::Zero(num_doc, num_topics);
    total_words_corpus[corpus] = 0.0;
  }

  for (int doc_id = 0; doc_id < num_doc; doc_id++) {
    corpus_ = corpus_id[doc_id];
    doc_w = W[doc_id];
    doc_len_ = doc_each_len[doc_id];
    for (int w_position = 0; w_position < doc_len_; w_position++) {
      w_ = doc_w[w_position];
      vocab_weights_corpus(corpus_, w_)++;
      (total_words_corpus[corpus_])++;
    }
  }

  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    if (weights_type == "inv-freq" || weights_type == "inv-freq-normalized") {
      vocab_weights_corpus.row(corpus) = (double)total_words_corpus[corpus] / (vocab_weights_corpus.row(corpus)).array();
    } else if (weights_type == "information-theory" || weights_type == "information-theory-normalized") {
      vocab_weights_corpus.row(corpus) = (vocab_weights_corpus.row(corpus)).array() / (double)total_words_corpus[corpus];
      vocab_weights_corpus.row(corpus) = (vocab_weights_corpus.row(corpus)).array().log();
      vocab_weights_corpus.row(corpus) = -(vocab_weights_corpus.row(corpus)).array() / log(2);
    }

    if (use_weights == 0) {
      Rcpp::Rcout << "Not using weights!! Check `options$use_weights`." << std::endl;
      vocab_weights_corpus.row(corpus) = RowVectorXd::Constant(num_vocab, 1.0);
    }
  }

  if (weights_type == "inv-freq-normalized" || weights_type == "information-theory-normalized") {
    NumericVector total_weights(num_corpus);
    int doc_len;
    for (int doc_id = 0; doc_id < num_doc; ++doc_id) {
      corpus_ = corpus_id[doc_id];
      doc_w = W[doc_id];
      doc_len = doc_each_len[doc_id];
      for (int w_position = 0; w_position < doc_len; ++w_position) {
        w_ = doc_w[w_position];
        total_weights[corpus_] += vocab_weights_corpus(corpus_, w_);
      }  
      for (int corpus = 0; corpus < num_corpus; ++corpus) {
        vocab_weights_corpus.row(corpus) = (vocab_weights_corpus.row(corpus)).array() * (double)total_words_corpus[corpus_] / total_weights[corpus_];
      }
    }
  }

  //
  // Construct (corpus specific) data matrices
  //
  vector<Triplet> trip_s1;  // for a sparse matrix
  total_words_weighted = 0.0;
  double temp;
  
  for (int doc_id = 0; doc_id < num_doc; ++doc_id) {
    corpus_ = corpus_id[doc_id];
    doc_s = S[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id], doc_c = C[doc_id];
    doc_len_ = doc_each_len[doc_id];
    
    for (int w_position = 0; w_position < doc_len_; ++w_position) {
      s_ = doc_s[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position], c_ = doc_c[w_position];
      if (s_ == 0){
        if (c_ == 0) {
          n_s0_c0_kv_all[corpus_](z_, w_) += vocab_weights_corpus(corpus_, w_);
          n_c0_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_s0_c0_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
          n_c0_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
        } else {
          n_c1_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_s0_c1_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_s0_c1_kv(z_, w_) += vocab_weights_corpus(corpus_, w_);
          n_c1_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
        }
        n_s0_kv(z_, w_) += vocab_weights_corpus(corpus_, w_);
        n_s0_k(z_) += vocab_weights_corpus(corpus_, w_);
      } else {
        trip_s1.push_back(Triplet(z_, w_, vocab_weights_corpus(corpus_, w_)));
        n_s1_k(z_) += vocab_weights_corpus(corpus_, w_);
        if (c_ == 0) {
          n_c0_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_c0_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
        } else {
          n_c1_k(z_) += vocab_weights_corpus(corpus_, w_);
          n_c1_k_all(z_, corpus_) += vocab_weights_corpus(corpus_, w_);
        }
      }
      n_dk(doc_id, z_) += vocab_weights_corpus(corpus_, w_);
      n_dk_all[corpus_](doc_id, z_) += vocab_weights_corpus(corpus_, w_);
    }
    
    temp = n_dk.row(doc_id).sum();
    doc_each_len_weighted.push_back(temp);
    total_words_weighted += temp;
  }
  n_s1_kv_multi.setFromTriplets(trip_s1.begin(), trip_s1.end());
}

void keyATMmulticov::initialize_specific()
{
  alpha = VectorXd::Zero(num_topics);  // used in iteration
  Alpha_all.resize(num_corpus);
  Lambda_all.resize(num_corpus);

  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    Alpha_all[corpus] = MatrixXd::Zero(num_doc_all[corpus], num_topics);
    Lambda_all[corpus] = MatrixXd::Zero(num_topics, num_cov_all[corpus]);
    for (int k = 0; k < num_topics; ++k) {
      for (int i = 0; i < num_cov_all[corpus]; ++i) {
        Lambda_all[corpus](k, i) = R::rnorm(mu_all[corpus], sigma_all[corpus]);
      }
    }
  }

  create_sufficient_stats();
}

void keyATMmulticov::resume_initialize_specific()
{
  // Alpha
  alpha = VectorXd::Zero(num_topics); // used in iteration
  Alpha_all.resize(num_corpus);
  Lambda_all.resize(num_corpus);

  // Lambda
  List Lambda_iter_all = stored_values["Lambda_iter"];
  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    Lambda_all[corpus] = MatrixXd::Zero(num_topics, num_cov_all[corpus]);
    List Lambda_iter = Lambda_iter_all[corpus];
    NumericMatrix Lambda_r = Lambda_iter[Lambda_iter.size() - 1];
    Lambda_all[corpus] = Rcpp::as<Eigen::MatrixXd>(Lambda_r);
  }

  create_sufficient_stats();
}

void keyATMmulticov::iteration_single(int it)
{ // Single iteration
  int doc_id_;
  int local_doc_id_;
  int corpus;
  int doc_length;
  int w_, z_, s_, c_;
  int new_z, new_s, new_c;
  int w_position;

  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle

  // Create Alpha for this iteration
  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    Alpha_all[corpus] = (X[corpus] * Lambda_all[corpus].transpose()).array().exp();
  }

  for (int ii = 0; ii < num_doc; ++ii) {
    doc_id_ = doc_indexes[ii];
    doc_s = S[doc_id_], doc_z = Z[doc_id_], doc_w = W[doc_id_], doc_c = C[doc_id_];
    corpus = corpus_id[doc_id_];
    doc_length = doc_each_len[doc_id_];
    local_doc_id_ = doc_lookup_global[doc_id_];
    
    alpha = Alpha_all[corpus].row(local_doc_id_).transpose(); // take out alpha

    token_indexes = sampler::shuffled_indexes(doc_length); //shuffle

    // Iterate each word in the document
    for (int jj = 0; jj < doc_length; ++jj) {
      w_position = token_indexes[jj];
      s_ = doc_s[w_position], z_ = doc_z[w_position], w_ = doc_w[w_position], c_ = doc_c[w_position];

      new_z = sample_z(alpha, z_, s_, w_, c_, doc_id_);
      doc_z[w_position] = new_z;

      if (keywords[new_z].find(w_) == keywords[new_z].end())
        continue;

      z_ = doc_z[w_position]; // use updated z
      new_s = sample_s(z_, s_, w_, c_, doc_id_);
      doc_s[w_position] = new_s;

      s_ = doc_s[w_position]; // use updated s
      new_c = sample_c(z_, s_, w_, c_, doc_id_);
      doc_c[w_position] = new_c;
    }

    Z[doc_id_] = doc_z;
    S[doc_id_] = doc_s;
    C[doc_id_] = doc_c;
  }

  sample_parameters(it);
}

void keyATMmulticov::sample_parameters(int it)
{
  sample_lambda();

  // Store lambda
  int r_index = it + 1;
  if (r_index % thinning == 0 || r_index == 1 || r_index == iter) {
    for (int corpus = 0; corpus < num_corpus; ++corpus) {
      Rcpp::NumericMatrix Lambda_R = Rcpp::wrap(Lambda_all[corpus]);
      List Lambda_iter_all = stored_values["Lambda_iter"];
      List Lambda_iter = Lambda_iter_all[corpus];
      Lambda_iter.push_back(Lambda_R);
      Lambda_iter_all[corpus] = Lambda_iter;
      stored_values["Lambda_iter"] = Lambda_iter_all;
    }
  }
}

int keyATMmulticov::sample_z(VectorXd &alpha, int z, int s, int w, int c, int doc_id)
{
  int new_z;
  double numerator, denominator;
  double sum;

  int corpus_ = corpus_id[doc_id];
  Ref<VectorXd> n_s0_c0_k_a = n_s0_c0_k_all.col(corpus_);

  // remove data
  if (s == 0) {
    n_s0_c0_k_a(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_k(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) -= vocab_weights_corpus(corpus_, w);
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) -= vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s1_k(z) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    }
  } else {
    Rcerr << "Error at sample_z, remove" << std::endl;
  }

  n_dk(doc_id, z) -= vocab_weights_corpus(corpus_, w);
  n_dk_all[corpus_](doc_id, z) -= vocab_weights_corpus(corpus_, w);

  new_z = -1; // initialize
  if (s == 0 && c == 0) {
    for (int k = 0; k < num_topics; ++k) {

      numerator = (beta + n_s0_c0_kv_all[corpus_](k, w)) *
        (n_s0_k(k) + prior_gamma(k, 1)) *
        (n_c0_k(k) + prior_omega(k, 1)) *
        (n_dk_all[corpus_](doc_id, k) + alpha(k));

      denominator = (Vbeta + n_s0_c0_k_a(k)) *
        (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1)) *
        (n_c1_k(k) + prior_omega(k, 0) + n_c0_k(k) + prior_omega(k, 1));

      z_prob_vec(k) = numerator / denominator;
    }

    sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  } else if (s == 0 && c == 1) {
    for (int k = 0; k < num_topics; ++k) {

      numerator = (beta_c + n_s0_c1_kv(k, w)) *
        (n_s0_k(k) + prior_gamma(k, 1)) *
        (n_c1_k(k) + prior_omega(k, 0)) *
        (n_dk_all[corpus_](doc_id, k) + alpha(k));

      denominator = (Vbeta_c + n_s0_c1_k(k)) *
        (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1)) *
        (n_c1_k(k) + prior_omega(k, 0) + n_c0_k(k) + prior_omega(k, 1));

      z_prob_vec(k) = numerator / denominator;
    }

    sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample
    
  } else if (s == 1 && c == 0) {
    for (int k = 0; k < num_topics; ++k) {
      if (keywords[k].find(w) == keywords[k].end()) {
        z_prob_vec(k) = 0.0;
        continue;
      } else {
        numerator = (beta_s + n_s1_kv_multi.coeffRef(k, w)) *
          (n_s1_k(k) + prior_gamma(k, 0)) *
          (n_c0_k(k) + prior_omega(k, 1)) *
          (n_dk_all[corpus_](doc_id, k) + alpha(k));

        denominator = (Lbeta_sk(k) + n_s1_k(k)) *
          (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1)) *
          (n_c1_k(k) + prior_omega(k, 0) + n_c0_k(k) + prior_omega(k, 1));

        z_prob_vec(k) = numerator / denominator;
      }
    }

    sum = z_prob_vec.sum(); // normalize
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  } else {
    for (int k = 0; k < num_topics; ++k) {
      if (keywords[k].find(w) == keywords[k].end()) {
        z_prob_vec(k) = 0.0;
        continue;
      } else {
        numerator = (beta_s + n_s1_kv_multi.coeffRef(k, w)) *
          (n_s1_k(k) + prior_gamma(k, 0)) *
          (n_c1_k(k) + prior_omega(k, 0)) *
          (n_dk_all[corpus_](doc_id, k) + alpha(k));
        denominator = (Lbeta_sk(k) + n_s1_k(k) ) *
          (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1)) *
          (n_c1_k(k) + prior_omega(k, 0) + n_c0_k(k) + prior_omega(k, 1));

        z_prob_vec(k) = numerator / denominator;
      }
    }

    sum = z_prob_vec.sum();
    new_z = sampler::rcat_without_normalize(z_prob_vec, sum, num_topics); // take a sample

  }

  // add back data counts
  if (s == 0) {
    n_s0_c0_k_a(new_z) += vocab_weights_corpus(corpus_, w);
    n_s0_k(new_z) += vocab_weights_corpus(corpus_, w);
    n_s0_kv(new_z, w) += vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](new_z, w) += vocab_weights_corpus(corpus_, w);
      n_c0_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_c0_k_all(new_z, corpus_) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(new_z, w) += vocab_weights_corpus(corpus_, w);
      n_c1_k_all(new_z, corpus_) += vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(new_z, w) += vocab_weights_corpus(corpus_, w);
    n_s1_k(new_z) += vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_c0_k_all(new_z, corpus_) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(new_z) += vocab_weights_corpus(corpus_, w);
      n_c1_k_all(new_z, corpus_) += vocab_weights_corpus(corpus_, w);
    }
  } else {
    Rcerr << "Error at sample_z, add" << std::endl;
  }
  n_dk_all[corpus_](doc_id, new_z) += vocab_weights_corpus(corpus_, w);
  n_dk(doc_id, new_z) += vocab_weights_corpus(corpus_, w);

  return new_z;
}

int keyATMmulticov::sample_s(int z, int s, int w, int c, int doc_id)
{
  int new_s;
  double numerator, denominator;
  double s0_prob;
  double s1_prob;
  double sum;

  // Get corpus id and choose containers
  int corpus_ = corpus_id[doc_id];
  Ref<VectorXd> n_s0_c0_k_a = n_s0_c0_k_all.col(corpus_);

  // remove data
  if (s == 0) {
    n_s0_c0_k_a(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_k(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) -= vocab_weights_corpus(corpus_, w);
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) -= vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s1_k(z) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
    }
  } else {
    Rcerr << "Error at sample_z, remove" << std::endl;
  }

  // newprob_s1()
  numerator = (beta_s + n_s1_kv_multi.coeffRef(z, w)) *
    ( n_s1_k(z) + prior_gamma(z, 0) );
  denominator = (Lbeta_sk(z) + n_s1_k(z) );
  s1_prob = numerator / denominator;

  // newprob_s0()
  if (c == 0) {
    numerator = (beta + n_s0_c0_kv_all[corpus_](z, w)) *
        (n_s0_k(z) + prior_gamma(z, 1));

    denominator = (Vbeta + n_s0_c0_k_a(z) );
  } else {
    numerator = (beta_c + n_s0_c1_kv(z, w)) *
        (n_s0_k(z) + prior_gamma(z, 1));

    denominator = (Vbeta_c + n_s0_c1_k(z) );
  }
  s0_prob = numerator / denominator;
  
  // Normalize
  sum = s0_prob + s1_prob;
  s1_prob = s1_prob / sum;

  new_s = R::runif(0,1) <= s1_prob;  //new_s = Bern(s0_prob, s1_prob);

  // add back data counts
  if (new_s == 0) {
    n_s0_c0_k_a(z) += vocab_weights_corpus(corpus_, w);
    n_s0_k(z) += vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) += vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) += vocab_weights_corpus(corpus_, w);
      n_c0_k(z) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) += vocab_weights_corpus(corpus_, w);
    }
  } else if (new_s == 1) {
    n_s1_kv_multi.coeffRef(z, w) += vocab_weights_corpus(corpus_, w);
    n_s1_k(z) += vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(z) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) += vocab_weights_corpus(corpus_, w);
    }
  }

  return new_s;
}

int keyATMmulticov::sample_c(int z, int s, int w, int c, int doc_id)
{
  int new_c;
  double numerator, denominator;
  double c0_prob;
  double c1_prob;
  double sum;

  // Get corpus id and choose containers
  int corpus_ = corpus_id[doc_id];
  Ref<VectorXd> n_s0_c0_k_a = n_s0_c0_k_all.col(corpus_);

  // remove data
  if (s == 0) {
    n_s0_c0_k_a(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_k(z) -= vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) -= vocab_weights_corpus(corpus_, w);
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) -= vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(z, w) -= vocab_weights_corpus(corpus_, w);
    n_s1_k(z) -= vocab_weights_corpus(corpus_, w);
    if (c == 0) {
      n_c0_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) -= vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) -= vocab_weights_corpus(corpus_, w);
    }
  } else {
    Rcerr << "Error at sample_c, remove" << std::endl;
  }

  // newprob_c1()
  if (s == 0) {
    numerator = (beta_c + n_s0_c1_kv(z, w)) *
      ( n_c1_k(z) + prior_omega(z, 0) );
    denominator = (Vbeta_c + n_s0_c1_k(z) );
  } else {
    numerator = (beta_s + n_s1_kv_multi.coeffRef(z, w)) *
      ( n_c1_k(z) + prior_omega(z, 0) );
    denominator = (Lbeta_sk(z) + n_s1_k(z) );
  }
  c1_prob = numerator / denominator;

  // newprob_c0()
  if (s == 0) {
    numerator = (beta + n_s0_c0_kv_all[corpus_](z, w)) *
        (n_c0_k(z) + prior_omega(z, 1));

    denominator = (Vbeta + n_s0_c0_k_a(z) );
  } else {
    numerator = (beta_c + n_s0_c1_kv(z, w)) *
        (n_c0_k(z) + prior_omega(z, 1));

    denominator = (Vbeta_c + n_s0_c1_k(z) );
  }
  c0_prob = numerator / denominator;
  
  // Normalize
  sum = c0_prob + c1_prob;
  c1_prob = c1_prob / sum;

  new_c = R::runif(0,1) <= c1_prob;

  // add back data counts
  if (s == 0) {
    n_s0_c0_k_a(z) += vocab_weights_corpus(corpus_, w);
    n_s0_k(z) += vocab_weights_corpus(corpus_, w);
    n_s0_kv(z, w) += vocab_weights_corpus(corpus_, w);
    if (new_c == 0) {
      n_s0_c0_kv_all[corpus_](z, w) += vocab_weights_corpus(corpus_, w);
      n_c0_k(z) += vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_s0_c1_kv(z, w) += vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) += vocab_weights_corpus(corpus_, w);
    }
  } else if (s == 1) {
    n_s1_kv_multi.coeffRef(z, w) += vocab_weights_corpus(corpus_, w);
    n_s1_k(z) += vocab_weights_corpus(corpus_, w);
    if (new_c == 0) {
      n_c0_k(z) += vocab_weights_corpus(corpus_, w);
      n_c0_k_all(z, corpus_) += vocab_weights_corpus(corpus_, w);
    } else {
      n_c1_k(z) += vocab_weights_corpus(corpus_, w);
      n_c1_k_all(z, corpus_) += vocab_weights_corpus(corpus_, w);
    }
  }

  return new_c;
}

void keyATMmulticov::sample_lambda()
{
  double start, end, previous_p, new_p, newlikelihood, slice_, current_lambda, store_loglik, newlambdallk;

  corpus_ids = sampler::shuffled_indexes(num_corpus);
  topic_ids = sampler::shuffled_indexes(num_topics);
  int k, t;
  const double A = slice_A;

  newlambdallk = 0.0;

  for (int jj = 0; jj < num_corpus; ++jj) {
    int corpus = corpus_ids[jj];
    cov_ids = sampler::shuffled_indexes(num_cov_all[corpus]);

    for (int kk = 0; kk < num_topics; ++kk) {
      k = topic_ids[kk];

      for (int tt = 0; tt < num_cov_all[corpus]; ++tt) {
        t = cov_ids[tt];
        store_loglik = likelihood_lambda(corpus, k, t);

        start = val_min; // shrinked value
        end = val_max; // shrinked value

        current_lambda = Lambda_all[corpus](k,t);
        previous_p = shrink(current_lambda, A);
        slice_ = store_loglik - std::log(A * previous_p * (1.0 - previous_p))
                + log(unif_rand()); // <-- using R random uniform

        for (int shrink_time = 0; shrink_time < max_shrink_time; ++shrink_time) {
          new_p = sampler::slice_uniform(start, end); // <-- using R function above
          Lambda_all[corpus](k,t) = expand(new_p, A); // expand

          newlambdallk = likelihood_lambda(corpus, k, t);

          newlikelihood = newlambdallk - std::log(A * new_p * (1.0 - new_p));

          if (slice_ < newlikelihood) {
            break;
          } else if (abs(end - start) < 1e-9) {
            Rcerr << "Shrinked too much. Using a current value." << std::endl;
            Lambda_all[corpus](k,t) = current_lambda;
            break;
          } else if (previous_p < new_p) {
            end = new_p;
          } else if (new_p < previous_p) {
            start = new_p;
          } else {
            Rcpp::stop("Something goes wrong in sample_lambda().");
          }

        }
      }
    }
  }
}

int keyATMmulticov::lookup_global_doc(int corpus, int id) {
  List corpus_lookup = doc_lookup[corpus];
  return corpus_lookup[id];
}

double keyATMmulticov::likelihood_lambda(int corpus, int k, int t)
{
  double loglik = 0.0;
  MatrixXd Alpha = (X[corpus] * Lambda_all[corpus].transpose()).array().exp();
  alpha = VectorXd::Zero(num_topics);

  for (int d = 0; d < num_doc_all[corpus]; ++d) {
    int global_d = lookup_global_doc(corpus, d);

    alpha = Alpha.row(d).transpose(); // Doc alpha, column vector

    loglik += mylgamma(alpha.sum());
    loglik -= mylgamma( doc_each_len_weighted[global_d] + alpha.sum() );

    loglik -= mylgamma(alpha(k));
    loglik += mylgamma( n_dk(global_d, k) + alpha(k) );
  }

  // Prior
  loglik += -0.5 * log(2.0 * PI_V * std::pow(sigma_all[corpus], 2.0) );
  loglik -= ( std::pow( (Lambda_all[corpus](k,t) - mu_all[corpus]) , 2.0) / (2.0 * std::pow(sigma_all[corpus], 2.0)) );

  return loglik;
}

double keyATMmulticov::loglik_total()
{
  double loglik = 0.0;
  double fixed_part = 0.0;

  for (int k = 0; k < num_topics; ++k) {
    for (int v = 0; v < num_vocab; ++v) { // word
      loglik += mylgamma(beta_c + n_s0_c1_kv(k, v)) - mylgamma(beta_c);
      for (int j = 0; j < num_corpus; ++j) {
        loglik += mylgamma(beta + n_s0_c0_kv_all[j](k, v)) - mylgamma(beta);
      }
    }

    loglik += mylgamma(beta_c * (double)num_vocab) - mylgamma(beta_c * (double)num_vocab + n_s0_c1_k(k));
    for (int j = 0; j < num_corpus; ++j) {
      loglik += mylgamma(beta * (double)num_vocab + n_s0_c0_k_all(k, j)) - mylgamma(beta * (double)num_vocab + n_s0_c0_k_all(k, j));
    }

    if (k < keyword_k) {
      for (SparseMatrix<double,RowMajor>::InnerIterator it(n_s1_kv_multi, k); it; ++it) {
        loglik += mylgamma(beta_s + it.value()) - mylgamma(beta_s);
      }
      loglik += mylgamma(beta_s * (double)keywords_num[k]) - mylgamma(beta_s * (double)keywords_num[k] + n_s1_k(k));
    
      loglik += mylgamma(prior_gamma(k, 0) + prior_gamma(k, 1)) - mylgamma(prior_gamma(k, 0)) - mylgamma(prior_gamma(k, 1));

      loglik += mylgamma(n_s0_k(k) + prior_gamma(k, 1))
                - mylgamma(n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1))
                + mylgamma(n_s1_k(k) + prior_gamma(k, 0));
    }

    for (int j = 0; j < num_corpus; ++j) {
      loglik += mylgamma(prior_omega(k, 0) + prior_omega(k, 1)) - mylgamma(prior_omega(k, 0)) - mylgamma(prior_omega(k, 1));
      
      loglik += mylgamma(n_c0_k_all(k, j) + prior_omega(k, 1))
                - mylgamma(n_c1_k_all(k, j) + prior_omega(k, 0) + n_c0_k_all(k, j) + prior_omega(k, 1))
                + mylgamma(n_c1_k_all(k, j) + prior_omega(k, 0));
    }

    for (int corpus = 0; corpus < num_corpus; ++corpus) {
      for (int t = 0; t < num_cov_all[corpus]; ++t) {
        loglik += -0.5 * log(2.0 * PI_V * std::pow(sigma_all[corpus], 2.0));
        loglik -= ( std::pow( (Lambda_all[corpus](k,t) - mu_all[corpus]) , 2.0) / (2.0 * std::pow(sigma_all[corpus], 2.0)) );
      }
    }
  }

  for (int corpus = 0; corpus < num_corpus; ++corpus) {
    Alpha_all[corpus] = (X[corpus] * Lambda_all[corpus].transpose()).array().exp();
    for (int d = 0; d < num_doc_all[corpus]; ++d) {
      alpha = Alpha_all[corpus].row(d).transpose();
      int global_d = lookup_global_doc(corpus, d);
      loglik += mylgamma(alpha.sum()) - mylgamma(doc_each_len_weighted[global_d] + alpha.sum());

      for (int k = 0; k < num_topics; ++k) {
        loglik += mylgamma(n_dk(global_d,k) + alpha(k)) - mylgamma(alpha(k));
      }
    }
  }
  
  return loglik;
}