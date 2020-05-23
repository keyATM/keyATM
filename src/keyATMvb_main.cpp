#include "keyATMvb_main.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


keyATMvb::keyATMvb(List model_)
{
  // Constructor
  model = model_;
}


keyATMvb::~keyATMvb()
{
  // Destructor
}


void keyATMvb::fit()
{
  // Read data, initialize the model and fit the model
  read_data();
  initialize();
  iteration();
}


void keyATMvb::read_data()
{
  // `common` reads data required in all models
  // `specific` reads model specific data
  read_data_specific();
  read_data_common();

  // Read words
  read_data_words();
}


void keyATMvb::read_data_common()
{
  // Read data
  W = model["W"]; Z = model["Z"];
  vocab = model["vocab"];
  model_name = Rcpp::as<std::string>(model["model"]);
  stored_values = model["stored_values"];

  num_doc = W.size();
  num_vocab = vocab.size();

  // Options
  options_list = model["options"];
  use_weight = options_list["use_weights"];
  vb_options = model["vb_options"];

  // Make an alpha matrix
  read_data_common_alpha();

  // Read beta
  beta = priors_list["beta"];
  Vbeta = num_vocab * beta;
  
  // Read gamma
  prior_gamma = MatrixXd::Zero(num_topics, 2);
  NumericMatrix RMatrix = priors_list["gamma"];
  prior_gamma = Rcpp::as<Eigen::MatrixXd>(RMatrix);

  // Weights
  NumericVector vocab_weights_R = stored_values["vocab_weights"];
  vocab_weights = Rcpp::as<Eigen::VectorXd>(vocab_weights_R);
}


void keyATMvb::read_data_specific()
{
  // For keyATM
  S = model["S"];
  regular_k = model["no_keyword_topics"];
  keywords_list = model["keywords"];
  keyword_k = keywords_list.size();
  num_topics = keyword_k + regular_k;

  priors_list = model["priors"];
  beta_s = priors_list["beta_s"];
}


void keyATMvb::read_data_common_alpha()
{
  // Make D \times K alpha matrix
  alphas = MatrixXd::Zero(num_doc, num_topics);

  if (model_name == "base") {
    read_data_common_alpha_base();
  } else if (model_name == "covariates") {
    read_data_common_alpha_cov();
  } else if (model_name == "dynamic") {
    read_data_common_alpha_hmm();
  } else {
    Rcpp::stop("Invalid model type");
  }

  alpha_d = alphas.rowwise().sum();
}


void keyATMvb::read_data_common_alpha_base()
{
  //
  // Reorganize alpha in the base model
  //
  List alpha_iter = stored_values["alpha_iter"];
  int total = alpha_iter.size();
  NumericVector alpha_R;
  VectorXd alpha = VectorXd::Zero(num_topics);

  int divide = 0;
  for (int iter = floor(total * 0.9); iter < total; iter++) {
    alpha_R = alpha_iter[iter];

    for (int k = 0; k < num_topics; k++) {
      alpha(k) += alpha_R[k]; 
    }

    divide += 1;
  }
  alpha = alpha.array() / (double)divide;  // take a mean

  // Into a D \times K form
  for (int d = 0; d < num_doc; d++) {
    alphas.row(d) = alpha.transpose(); 
  }
}


void keyATMvb::read_data_common_alpha_cov()
{
  //
  // Reorganize alpha in the cov model
  //
  List model_settings = model["model_settings"];
  NumericMatrix used_data_R = model_settings["covariates_data_use"];
  MatrixXd use_data = Rcpp::as<Eigen::MatrixXd>(used_data_R);
  int num_cov = use_data.cols();

  List Lambda_iter = stored_values["Lambda_iter"];
  int total = Lambda_iter.size();
  MatrixXd Lambda = MatrixXd::Zero(num_topics, num_cov);

  int divide = 0;
  for (int iter = floor(total * 0.9); iter < total; iter++) {
    NumericMatrix Lambda_R = Lambda_iter[iter];

    for (int k = 0; k < num_topics; k++) {
      for (int c = 0; c < num_cov; c++) {
        Lambda(k, c) += Lambda_R(k, c); 
      } 
    }

    divide += 1;
  }
  Lambda = Lambda.array() / (double)divide;  // take a mean

  // Get Alpha
  alphas = (use_data * Lambda.transpose()).array().exp();
}


void keyATMvb::read_data_common_alpha_hmm()
{
  //
  // Reorganize alpha in the HMM model
  //

  // R
  List R_iter = stored_values["R_iter"];
  int last = R_iter.size() - 1;
  NumericVector R_R = R_iter[last];
  VectorXd R = Rcpp::as<Eigen::VectorXd>(R_R);

  // Time index
  List model_settings = model["model_settings"];
  NumericVector time_index_R = model_settings["time_index"];
  VectorXd time_index = Rcpp::as<Eigen::VectorXd>(time_index_R);
  time_index = time_index.array() - 1;

  // alpha
  List alpha_iter = stored_values["alpha_iter"];
  int total = alpha_iter.size();
  int num_states = R.maxCoeff() + 1;  // adjust index
  MatrixXd alpha_t = MatrixXd::Zero(num_states, num_topics);

  int divide = 0;
  for (int iter = floor(total * 0.9); iter < total; iter++) {
    NumericMatrix alpha_R = alpha_iter[iter];

    for (int s = 0; s < num_states; s++) {
      for (int k = 0; k < num_topics; k++) {
        alpha_t(s, k) += alpha_R(s, k);
      }
    }

    divide += 1;
  }
  alpha_t = alpha_t.array() / (double)divide;

  // Into a D \times K form
  int time, state;
  for (int d = 0; d < num_doc; d++) {
    time =  time_index(d);
    state = R(time);
    alphas.row(d) = alpha_t.row(state);
  }

}


void keyATMvb::read_data_words()
{
  read_data_keywords();  // for keyATM

}


void keyATMvb::read_data_keywords()
{
  int wd_id;
  IntegerVector wd_ids;
  for (int ii = 0; ii < keyword_k; ii++) {
    wd_ids = keywords_list[ii];
    keywords_num.push_back(wd_ids.size());
    
    std::unordered_set<int> keywords_set;
    for (int jj = 0; jj < wd_ids.size(); jj++) {
      wd_id = wd_ids(jj);
      keywords_set.insert(wd_id);
      keywords_all.insert(wd_id);
    }
    
    keywords.push_back(keywords_set);
  }

  for (int i = keyword_k; i < num_topics; i++) {
    std::unordered_set<int> keywords_set{ -1 };
  
    keywords_num.push_back(0);
    keywords.push_back(keywords_set);
  }
}




void keyATMvb::initialize()
{
  initialize_common();
  initialize_specific();
}


void keyATMvb::initialize_common()
{
  initialize_weightedlen();
  initialize_common_MCMCcount();
  initialize_common_q();
  initialize_common_expectation();

  // During the iteration
  z_prob_vec = VectorXd::Zero(num_topics);
  s_prob_vec = VectorXd::Zero(2);
  s0_temp = VectorXd::Zero(num_topics);
  s1_temp = VectorXd::Zero(num_topics);


  // Initialize weights
}


void keyATMvb::initialize_common_MCMCcount()
{
  // Counts from MCMC used in initialize_common_q
  int w, z, s;
  int doc_len;

  n_s0_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_s1_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_s0_k = VectorXd::Zero(num_topics);
  n_s1_k = VectorXd::Zero(num_topics);
  n_dk = MatrixXd::Zero(num_doc, num_topics);
  
  for (int doc_id = 0; doc_id < num_doc; doc_id++) {
    doc_s = S[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
    doc_len = doc_each_len[doc_id];

    for (int w_position = 0; w_position < doc_len; w_position++) {
      s = doc_s[w_position], z = doc_z[w_position], w = doc_w[w_position];
      if (s == 0){
        n_s0_kv(z, w) += vocab_weights(w);
        n_s0_k(z) += vocab_weights(w);
      } else {
        n_s1_kv(z, w) += vocab_weights(w);
        n_s1_k(z) += vocab_weights(w);
      }
      n_dk(doc_id, z) += 1.0;
    }
  }
  
}


void keyATMvb::initialize_common_q()
{
  // Initialization using MCMC results

  List doc_z;
  List doc_s;
  int s, z, w;
  int doc_len;
  
  // Decide the values
  z_prob_vec = VectorXd::Zero(num_topics);

  for (int doc_id = 0; doc_id < num_doc; doc_id++) {
    doc_len = doc_each_len[doc_id];
    std::vector<std::vector<double>> qzd;  // store document level qz
    doc_z = Z[doc_id];

    std::vector<std::vector<double>> qsd;
    doc_s = S[doc_id];
    doc_w = W[doc_id];

    for (int w_position = 0; w_position < doc_len; w_position++) {
      z = doc_z[w_position];
      s = doc_s[w_position];
      w = doc_w[w_position];

      // Initialize qz
      std::vector<double> qzdk(num_topics, 0.0);
      initialize_common_qz(doc_id, w, z, s, qzdk);
      qzd.push_back(qzdk);
      
      // Initialize qs
      std::vector<double> qsds(2, 0.0);
      initialize_common_qs(doc_id, w, z, s, qsds);
      qsd.push_back(qsds);

    }

    qz.push_back(qzd);
    qs.push_back(qsd);
  }
}


void keyATMvb::initialize_common_qz(int doc_id, int w, int z, int s, vector<double> &qzdk)
{
  double numerator, denominator;
  double sum;

  if (s == 0) {
    for (int k = 0; k < num_topics; ++k) {

      numerator = (beta + n_s0_kv(k, w)) *
        (n_s0_k(k) + prior_gamma(k, 1)) *
        (n_dk(doc_id, k) + alphas(doc_id, k));

      denominator = ((double)num_vocab * beta + n_s0_k(k)) *
        (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1));

      z_prob_vec(k) = numerator / denominator;
      // if (k == z)  // checking, debug
      //   z_prob_vec(k) *= 5;
    }

  } else {
    for (int k = 0; k < num_topics; ++k) {
      if (keywords[k].find(w) == keywords[k].end()) {
        z_prob_vec(k) = 0.0;
        continue;
      } else { 
        numerator = (beta_s + n_s1_kv.coeffRef(k, w)) *
          (n_s1_k(k) + prior_gamma(k, 0)) *
          (n_dk(doc_id, k) + alphas(doc_id, k));
      }
      denominator = ((double)keywords_num[k] * beta_s + n_s1_k(k) ) *
        (n_s1_k(k) + prior_gamma(k, 0) + n_s0_k(k) + prior_gamma(k, 1));

      z_prob_vec(k) = numerator / denominator;
      // if (k == z)  // checking, debug
      //   z_prob_vec(k) *= 5;
    }

  }

  sum = z_prob_vec.sum();
  z_prob_vec = z_prob_vec.array() / sum;

  for (int k = 0; k < num_topics; k++) {
    qzdk[k]  = z_prob_vec(k);
  }

}


void keyATMvb::initialize_common_qs(int doc_id, int w, int z, int s, vector<double> &qsds)
{
  double numerator, denominator;
  double s0_prob, s1_prob;
  double sum;

  if (z >= keyword_k) {
    // Not a keyword topic
    qsds[0] = 0.99999;
    qsds[1] = 0.00001;
  } else if (keywords[z].find(w) == keywords[z].end()) {
    // Not a keyword
    qsds[0] = 0.99999;
    qsds[1] = 0.00001;
  } else {
    // One of the keywords
    numerator = (beta_s + n_s1_kv.coeffRef(z, w)) *
      ( n_s1_k(z) + prior_gamma(z, 0) );
    denominator = ((double)keywords_num[z] * beta_s + n_s1_k(z) );
    s1_prob = numerator / denominator;

    // newprob_s0()
    numerator = (beta + n_s0_kv(z, w)) *
      (n_s0_k(z) + prior_gamma(z, 1));

    denominator = ((double)num_vocab * beta + n_s0_k(z) );
    s0_prob = numerator / denominator;

    // Normalize
    sum = s0_prob + s1_prob;
    s0_prob = s0_prob / sum;
    s1_prob = s1_prob / sum;

    qsds[0] = s0_prob;
    qsds[1] = s1_prob;
  }
}


void keyATMvb::initialize_common_expectation()
{
  // Calculate a lot of expectations
  int w, z, s;
  int doc_len;

  n_s0_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_s1_kv = MatrixXd::Zero(num_topics, num_vocab);
  n_s0_k = VectorXd::Zero(num_topics);
  n_s1_k = VectorXd::Zero(num_topics);
  n_dk = MatrixXd::Zero(num_doc, num_topics);

  
  for (int doc_id = 0; doc_id < num_doc; doc_id++) {
    doc_s = S[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
    doc_len = doc_each_len[doc_id];
  
    for (int w_position = 0; w_position < doc_len; w_position++) {
      s = doc_s[w_position], z = doc_z[w_position], w = doc_w[w_position];
    
      for (int k = 0; k < num_topics; k++){
        n_s0_kv(k, w) += qz[doc_id][w_position][k] * qs[doc_id][w_position][0] * vocab_weights(w);
        n_s1_kv(k, w) += qz[doc_id][w_position][k] * qs[doc_id][w_position][1] * vocab_weights(w);

        n_s0_k(k) += qz[doc_id][w_position][k] * qs[doc_id][w_position][0] * vocab_weights(w);
        n_s1_k(k) += qz[doc_id][w_position][k] * qs[doc_id][w_position][1] * vocab_weights(w);

        n_dk(doc_id, k) += qz[doc_id][w_position][k];
      }
    }
  }
}


void keyATMvb::initialize_weightedlen()
{
  // Construct vocab weights
  int doc_len;
  double doc_len_weighted;
  total_words = 0;
  total_words_weighted = 0.0;
  IntegerVector doc_s, doc_z, doc_w;
  int w;
  
  for (int doc_id = 0; doc_id < num_doc; doc_id++) {
    doc_w = W[doc_id];
    doc_len = doc_w.size();
    doc_each_len.push_back(doc_len);
  
    doc_len_weighted = 0.0;
    for (int w_position = 0; w_position < doc_len; w_position++) {
      w = doc_w[w_position];
      total_words_weighted += vocab_weights(w);
      doc_len_weighted += vocab_weights(w);
    }
    doc_each_len_weighted.push_back(doc_len_weighted);

    total_words += doc_len;
  }

}


void keyATMvb::initialize_specific()
{

}


void keyATMvb::iteration()
{
  // Convergence settings
  double convtol = vb_options["convtol"];
  double change_rate = 1.00;
  int max_VB_iter = 5000;
  int doc_id;

  // Perplexity settings
  double perplexity;
  double previous_perplexity = -100.0;
  num_doc_perp = min((int)ceil(num_doc * 0.1), 100);
  ppl_doc_indexes = sampler::shuffled_indexes(num_doc_perp);

  if (num_doc_perp == num_doc) {
    ppl_words = total_words; 
  } else {
    // If you check a subset of documents
    ppl_words = 0.0;
    for (int d = 0; d < num_doc_perp; d++) {
      doc_id = ppl_doc_indexes[d];
      ppl_words += doc_each_len[doc_id];
    }
  }

  //
  // Iteration
  //
  int count = 1;
  while (change_rate > convtol) {
    iteration_single(); 
  
    if (count % 1 == 0) {
      perplexity = calc_perplexity(count);

      if (previous_perplexity < 0) {
        previous_perplexity = perplexity; 
      } else {
        change_rate = (previous_perplexity - perplexity) / previous_perplexity;
        previous_perplexity = perplexity;
      }
      Rcpp::Rcout << "Perplexity [" << count << "]: " << perplexity << " / ";
      Rcpp::Rcout << "Convergence [" << count << "]: " << change_rate << endl;
    }
  
  
    // Check keybord interruption to cancel the iteration
    checkUserInterrupt();

    count += 1;
    if (count == max_VB_iter)
      break;
  }
}


void keyATMvb::iteration_single()
{
  update_q();
}


void keyATMvb::update_q()
{
  doc_indexes = sampler::shuffled_indexes(num_doc); // shuffle
  int doc_id;
  int doc_len;
  double temp_sum;
  int v;

  for (int ii = 0; ii < num_doc; ii++) {
    doc_id = doc_indexes[ii];
    doc_w = W[doc_id];
    doc_len = doc_each_len[doc_id];
  
    for (int w_position = 0; w_position < doc_len; w_position++) {
      v = doc_w[w_position];
      update_decrese_count(doc_id, w_position, v);
    
      //
      // Update qz
      //
      for (int k = 0; k < num_topics; k++) {
        z_prob_vec(k) = exp(
                          qs[doc_id][w_position][0] * 
                              ( 
                               // log(n_s0_kv(k, v) + beta) - log(n_s0_k(k) + Vbeta)
                               //   + log(n_s0_k(k) + prior_gamma(k, 1))
                               //   - log(n_s0_k(k) + n_s1_k(k) + prior_gamma(k, 0) + prior_gamma(k, 1))
                               log( (n_s0_kv(k, v) + beta) / (n_s0_k(k) + Vbeta)
                                     * (n_s0_k(k) + prior_gamma(k, 1))
                                     / (n_s0_k(k) + n_s1_k(k) + prior_gamma(k, 0) + prior_gamma(k, 1)) )
                              )
    
                          + qs[doc_id][w_position][1] * 
                              ( 
                               // log(n_s1_kv(k, v) + beta_s) - log(n_s1_k(k) + (double)keywords_num[k] * beta_s)
                               //   + log(n_s1_k(k) + prior_gamma(k, 0))
                               //   - log(n_s0_k(k) + n_s1_k(k) + prior_gamma(k, 0) + prior_gamma(k, 1))
                               log( (n_s1_kv(k, v) + beta_s) / (n_s1_k(k) + (double)keywords_num[k] * beta_s)
                                     * (n_s1_k(k) + prior_gamma(k, 0))
                                     / (n_s0_k(k) + n_s1_k(k) + prior_gamma(k, 0) + prior_gamma(k, 1)) )
                              )
                        )
                        * (n_dk(doc_id, k) + alphas(doc_id, k));
      }
      
      temp_sum = z_prob_vec.sum();
      for (int k = 0; k < num_topics; k++) {
        // Normalize
        qz[doc_id][w_position][k] = z_prob_vec(k) / temp_sum; 
      }

 
      //
      // Update qs
      //
      if (keywords_all.find(v) == keywords_all.end()) {
        // No chance of using keywords distribution
        qs[doc_id][w_position][0] = 0.99999;
        qs[doc_id][w_position][1] = 0.00001;
      } else {
        // Can use a keyword
      
        for (int k = 0; k < num_topics; k++) {
          s0_temp(k) = qz[doc_id][w_position][k] * 
                          (
                            // log(n_s0_kv(k, v) + beta) - log(n_s0_k(k) + Vbeta)
                            //   + log(n_s0_k(k) + prior_gamma(k, 1))
                            log( (n_s0_kv(k, v) + beta) / (n_s0_k(k) + Vbeta)
                                  * (n_s0_k(k) + prior_gamma(k, 1)) )
                          );
      
          s1_temp(k) = qz[doc_id][w_position][k] * 
                          (
                            // log(n_s1_kv(k, v) + beta_s) - log(n_s1_k(k) + (double)keywords_num[k] * beta_s)
                              // + log(n_s1_k(k) + prior_gamma(k, 0))
                            log( (n_s1_kv(k, v) + beta_s) / (n_s1_k(k) + (double)keywords_num[k] * beta_s)
                                  *  (n_s1_k(k) + prior_gamma(k, 0)) )
                          );
                     
        }
      
        s_prob_vec(0) = exp(s0_temp.sum());
        s_prob_vec(1) = exp(s1_temp.sum());
      
        temp_sum = s_prob_vec.sum();
        for (int s = 0; s < 2; s++) {
          // Normalize 
          qs[doc_id][w_position][s] = s_prob_vec(s) / temp_sum;
        }
      
      }

      update_increase_count(doc_id, w_position, v);
    }

  }

}


void keyATMvb::update_decrese_count(int doc_id, int w_position, int v)
{
  double temp0;
  double temp1;

  for(int k = 0; k < num_topics; k++) {
    temp0 = qz[doc_id][w_position][k] * qs[doc_id][w_position][0] * vocab_weights(v);
    temp1 = qz[doc_id][w_position][k] * qs[doc_id][w_position][1] * vocab_weights(v);
    n_s0_kv(k, v) -= temp0;
    n_s1_kv(k, v) -= temp1;

    n_s0_k(k) -= temp0;
    n_s1_k(k) -= temp1;

    n_dk(doc_id, k) -= qz[doc_id][w_position][k];
  }
}


void keyATMvb::update_increase_count(int doc_id, int w_position, int v)
{
  double temp0;
  double temp1;

  for(int k = 0; k < num_topics; k++) {
    temp0 = qz[doc_id][w_position][k] * qs[doc_id][w_position][0] * vocab_weights(v);
    temp1 = qz[doc_id][w_position][k] * qs[doc_id][w_position][1] * vocab_weights(v);

    n_s0_kv(k, v) += temp0;
    n_s1_kv(k, v) += temp1;

    n_s0_k(k) += temp0;
    n_s1_k(k) += temp1;

    n_dk(doc_id, k) += qz[doc_id][w_position][k];
  }
}


double keyATMvb::calc_perplexity(int iter)
{
  double llk = 0.0;
  double llk_v = 0.0;
  double ppl;
  int doc_id;
  int doc_len;
  int v;

  for (int ii = 0; ii < num_doc_perp; ii++) {
    // Limit the number of documents to calculate the perplexity
    // to reduce the computation time

    doc_id = ppl_doc_indexes[ii];
    doc_len = doc_each_len[doc_id];
    doc_w = W[doc_id];

    for (int i = 0; i < doc_len; i++) {
      v = doc_w[i];
      llk_v = 0.0;

      // Word
      for (int k = 0; k < num_topics; k++) {
        llk_v += (
                   (n_s0_kv(k, v) + beta) / (n_s0_k(k) + Vbeta) 
                    * (n_s0_k(k) + prior_gamma(k, 1)) / (n_s0_k(k) + prior_gamma(k, 1) + n_s1_k(k) + prior_gamma(k, 0))
                  +
                   (n_s1_kv(k, v) + beta_s) / (n_s1_k(k) + (double)keywords_num[k] * beta_s)
                    * (n_s1_k(k) + prior_gamma(k, 0)) / (n_s0_k(k) + prior_gamma(k, 1) + n_s1_k(k) + prior_gamma(k, 0))
                 )
                  * (n_dk(doc_id, k) + alphas(doc_id, k)) / (doc_each_len_weighted[doc_id] + alpha_d(doc_id));
      }
      llk += log(llk_v);
    }
  
  }

  ppl = exp(-llk / ppl_words);

  // Store perplexity
  store_perplexity(iter, ppl);

  return ppl;
}


void keyATMvb::store_perplexity(int iter, double ppl)
{
  Perplexity_value.push_back(ppl);
  Perplexity_iter.push_back(iter);
  Perplexity["value"] = Perplexity_value;
  Perplexity["iter"] = Perplexity_iter;

  vb_options["Perplexity_VB"] = Perplexity;
}


List keyATMvb::return_model()
{
  // Return output to R
  get_QOI();

  // Organize values
  model["options"] = options_list;
  model["vb_options"] = vb_options;
  return model;
}

void keyATMvb::get_QOI()
{
  // Get z and s from qz and qs, respectively, and use it as the last iteration
  int z, s;
  int doc_len;

  for (int doc_id = 0; doc_id < num_doc; doc_id++) {
    doc_len = doc_each_len[doc_id];
    doc_z = Z[doc_id];
    doc_s = S[doc_id];

    for (int i = 0; i < doc_len; i++){
      z = max_element(qz[doc_id][i].begin(), qz[doc_id][i].end()) - qz[doc_id][i].begin();
      s = max_element(qs[doc_id][i].begin(), qs[doc_id][i].end()) - qs[doc_id][i].begin();
    
      doc_z[i] = z;
      doc_s[i] = s;
    }

    Z[doc_id] = doc_z;
    S[doc_id] = doc_s;
  }
}

