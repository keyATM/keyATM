#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;

double logsumexp (double &x, double &y, bool flg){
  if (flg) return y; // init mode
  if (x == y) return x + 0.69314718055; // log(2)
  double vmin = std::min (x, y);
  double vmax = std::max (x, y);
  if (vmax > vmin + 50) {
    return vmax;
  } else {
    return vmax + std::log (std::exp (vmin - vmax) + 1.0);
  }
}

double logsumexp_Eigen(Eigen::VectorXd &vec){
  double sum = 0.0;
  int index;
  for(size_t i = 0; i < vec.size(); i++){
    index = static_cast<int>(i);
    sum = logsumexp(sum, vec[index], (index == 0));
  }
  return sum;
}

double gammapdfln(double x, double a, double b){
 return a * log(b) - lgamma(a) + (a-1.0) * log(x) - b * x;
}

int rcat(Eigen::VectorXd &prob){
  // Multi(x, 1), return category index
  double u = R::runif(0, 1);
  double temp = 0.0;
  int index = 0;
  for(int ii = 0; ii < prob.size(); ii++){
    temp += prob(ii);
    if(u < temp){
      index = ii;
      break;
    }
  }
  return index;
}

// int bern(double &prob0, double &prob1){
//  double value = R::runif(0,1);
//  if(value <= prob1)
//    return 1;
//  return 0;
//}
//
// This function replaced by:   R::runif(0,1) <= prob1

int sample_z(Eigen::SparseMatrix<int, Eigen::RowMajor>& n_x0_kv,
             Eigen::SparseMatrix<int, Eigen::RowMajor>& n_x1_kv,
             Eigen::VectorXi& n_x0_k,
             Eigen::VectorXi& n_x1_k,
             Eigen::MatrixXd& n_dk,
             Eigen::VectorXd& alpha,
             std::vector<int>& seed_num,
             int x, int z, int w, int doc_id, int num_vocab,
             int num_topics, int k_seeded,
             double gamma_1, double gamma_2, double beta, double beta_s,
             std::vector< std::unordered_map<int, double> > & phi_s){
  // remove data
  if (x == 0){
    n_x0_kv.coeffRef(z, w) -= 1;
    n_x0_k(z) -= 1;
  } else {
    n_x1_kv.coeffRef(z, w) -= 1;
    n_x1_k(z) -= 1;
  }
  n_dk.coeffRef(doc_id, z) -= 1;

  // draw z
  VectorXd z_prob_vec = VectorXd::Zero(num_topics);
  int new_z = -1; // debug
  double numerator, denominator;
  if (x == 0){
    for (int k = 0; k < num_topics; k++){
      numerator = log(beta + (double)n_x0_kv.coeffRef(k, w)) +
        log((double)n_x0_k(k) + gamma_2) +
        log((double)n_dk.coeffRef(doc_id, k) + alpha(k));
      denominator = log((double)num_vocab * beta + (double)n_x0_kv.row(k).sum()) +
        log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
      z_prob_vec(k) = numerator - denominator;
    }
    double sum = logsumexp_Eigen(z_prob_vec); // normalize
    for (int k = 0; k < num_topics; k++)
      z_prob_vec(k) = exp(z_prob_vec(k) - sum);
    new_z = rcat(z_prob_vec); // take a sample

  } else {
    std::vector<int> make_zero_later;
    for (int k = 0; k < k_seeded; k++){ // <--------- NOTE: using k_seeded
      if (phi_s[k].find(w) == phi_s[k].end()){
        z_prob_vec(k) = 1.0;
        make_zero_later.push_back(k);
        continue;
      } else{ // w not one of the seeds
        numerator = log(beta_s + (double)n_x1_kv.coeffRef(k, w)) +
          log( ((double)n_x1_k(k) + gamma_1) ) +
          log( ((double)n_dk.coeffRef(doc_id, k) + alpha(k)) );
      }
      denominator = log((double)seed_num[k] * beta_s + (double)n_x1_kv.row(k).sum() ) +
        log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
      z_prob_vec(k) = numerator - denominator;
    }
    double sum = logsumexp_Eigen(z_prob_vec);
    for (int k = 0; k < num_topics; k++)
      z_prob_vec(k) = exp(z_prob_vec(k) - sum);
    for (int k = 0; k < make_zero_later.size(); k++) // zero out elements
      z_prob_vec(k) = 0.0;
    z_prob_vec = z_prob_vec / z_prob_vec.sum(); // and renormalize
    new_z = rcat(z_prob_vec); // take a sample
  }

  // add back data counts
  if (x == 0){
    n_x0_kv.coeffRef(z, w) += 1;
    n_x0_k(z) += 1;
  } else {
    n_x1_kv.coeffRef(z, w) += 1;
    n_x1_k(z) += 1;
  }
  n_dk.coeffRef(doc_id, z) += 1;

  return new_z;
}

int sample_x(SparseMatrix<int, RowMajor>& n_x0_kv,
             SparseMatrix<int, RowMajor>& n_x1_kv,
             VectorXi& n_x0_k,
             VectorXi& n_x1_k,
             MatrixXd& n_dk,
             VectorXd& alpha,
             std::vector<int>& seed_num,
             int x, int z, int w, int doc_id, int num_vocab,
             int num_topics, int k_seeded,
             double gamma_1, double gamma_2, double beta, double beta_s,
             std::vector< std::unordered_map<int, double> > & phi_s){
  // remove data
  if (x == 0){
    n_x0_kv.coeffRef(z, w) -= 1;
    n_x0_k(z) -= 1;
  } else {
    n_x1_kv.coeffRef(z, w) -= 1;
    n_x1_k(z) -= 1;
  }
  n_dk.coeffRef(doc_id, z) -= 1;

  // newprob_x1()
  double x1_logprob;
  int k = z;
  double numerator;
  double denominator;
  if(phi_s[k].find(w) == phi_s[k].end()){ // not a seed word
    x1_logprob = -1.0;
  } else {
    numerator = log(beta_s + (double)n_x1_kv.coeffRef(k, w)) +
      log( ((double)n_x1_k(k) + gamma_1) );
    denominator = log((double)seed_num[k] * beta_s + (double)n_x1_kv.row(k).sum() ) +
      log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
    x1_logprob = numerator - denominator;
  }

  int new_x;
  if(x1_logprob == -1.0){
    // if probability of x_di = 1 case is 0, it should be x=0 (regular topic)
    new_x = 0;
  } else {
    // newprob_x0()
    int k = z;
    numerator = log(beta + (double)n_x0_kv.coeffRef(k, w)) +
      log((double)n_x0_k(k) + gamma_2);
    denominator = log((double)num_vocab * beta + (double)n_x0_kv.row(k).sum() ) +
      log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
    double x0_logprob = numerator - denominator;

    // Normalize
    double sum = 0.0;
    double prob_array[] = {x0_logprob, x1_logprob};
    for (int i = 0; i < 2; i++)
      sum = logsumexp(sum, prob_array[i], (i == 0));

    double x0_prob = exp( x0_logprob - sum);
    double x1_prob = exp( x1_logprob - sum);
    new_x = R::runif(0,1) <= x1_prob;  //new_x = Bern(x0_prob, x1_prob);
  }

  // add back data counts
  if (x == 0){
    n_x0_kv.coeffRef(z, w) += 1;
    n_x0_k(z) += 1;
  } else {
    n_x1_kv.coeffRef(z, w) += 1;
    n_x1_k(z) += 1;
  }
  n_dk.coeffRef(doc_id, z) += 1;

  return new_x;
}


//' Run the Gibbs sampler
//'
//' @param model A model, from \code{init} or a previous invocation of \code{train}
//' @param k_seeded How many topics are seeded
//' @param k_free How many regular unseeded topics are required
//' @param alpha_k A starting value for alpha (will be removed when alpha updates are back in)
//' @param iter Required number of iterations
//'
//' @export
// [[Rcpp::export]]
List train(List model, int k_seeded, int k_free, double alpha_k, int iter = 0){

  List W = model["W"], Z = model["Z"], X = model["X"];
  StringVector files = model["files"], vocab = model["vocab"];
  List seeds = model["seeds"]; // Now convert this to T&S's phi_s format

  std::vector< std::unordered_map<int, double> > phi_s(seeds.size());
  std::vector<int> seed_num(seeds.size());
  for (int ii = 0; ii < seeds.size(); ii++){
    IntegerVector wd_ids = seeds[ii];
    seed_num[ii] = wd_ids.size();
    std::unordered_map<int, double> phi_sk;
    for (int jj = 0; jj < wd_ids.size(); jj++)
      phi_sk[wd_ids(jj)] = 1.0 / wd_ids.size();
    phi_s[ii] = phi_sk;
  }

  // alpha-related constants
  int num_topics = k_seeded + k_free;
  VectorXd alpha = VectorXd::Constant(num_topics, alpha_k);

  // document-related constants
  int num_vocab = vocab.size(), num_doc = files.size();

  // distributional constants
  double gamma_1 = 1.0, gamma_2 = 1.0;
  double lambda_1 = 1.0, lambda_2 = 2.0;
  double beta = 0.01, beta_s = 0.1;

  // storage for sufficient statistics and their margins
  SparseMatrix<int, RowMajor> n_x0_kv = SparseMatrix<int, RowMajor> (num_topics, num_vocab);
  SparseMatrix<int, RowMajor> n_x1_kv = SparseMatrix<int, RowMajor> (num_topics, num_vocab);
  MatrixXd n_dk = MatrixXd::Zero(num_doc, num_topics);
  MatrixXd theta_dk = MatrixXd::Zero(num_doc, num_topics); // Does this need initialization?
  VectorXi n_x0_k = VectorXi::Zero(num_topics);
  VectorXi n_x1_k = VectorXi::Zero(num_topics);

  for(int doc_id = 0; doc_id < num_doc; doc_id++){
    IntegerVector doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
    for(int w_position = 0; w_position < doc_x.size(); w_position++){
      int x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];
      if (x == 0){
        n_x0_kv.coeffRef(z, w) += 1;
        n_x0_k(z) += 1;
      } else {
        n_x1_kv.coeffRef(z, w) += 1;
        n_x1_k(z) += 1;
      }
      n_dk.coeffRef(doc_id, z) += 1.0;
    }
  }
  int total_words = n_dk.sum();

  // Sampler: No randomized update sequence this time
  for (int ii = 0; ii < iter; ii++){
    for (int doc_id = 0; doc_id < num_doc; doc_id++){
      IntegerVector doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
      for (int w_position = 0; w_position < doc_x.size(); w_position++){
        int x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];

        doc_z[w_position] = sample_z(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk,
                                     alpha,seed_num, x, z, w, doc_id,
                                     num_vocab, num_topics, k_seeded, gamma_1,
                                     gamma_2, beta, beta_s, phi_s);

        doc_x[w_position] = sample_z(n_x0_kv, n_x1_kv, n_x0_k, n_x1_k, n_dk,
                                    alpha,seed_num, x, z, w, doc_id,
                                    num_vocab, num_topics, k_seeded, gamma_1,
                                    gamma_2, beta, beta_s, phi_s);
      }
    }
    // update_alpha();

    // log likelihood
    double prod_k = 0.0;// calc_loglik_prod_k();
    for(int k=0; k<num_topics; ++k){
      // (a) Seed Topic Part
      prod_k += lgamma((double)num_vocab * beta_s); // first term numerator
      prod_k -= lgamma(beta_s) * (double)num_vocab; // first term denominator
      prod_k -= lgamma( (double)num_vocab * beta_s + (double)n_x1_kv.row(k).sum()  ); // second term denominator
      // Regular Topic Part
      prod_k += lgamma((double)num_vocab * beta); //probably constant
      prod_k -= lgamma(beta) * (double)num_vocab; // probably constant
      prod_k -= lgamma( (double)num_vocab * beta + (double)n_x0_kv.row(k).sum()  ); // last part denominator
      // (b) second part
      prod_k += lgamma((double)n_x1_k(k) + gamma_1) + lgamma((double)n_x0_k(k) + gamma_2);
      prod_k -=  lgamma( (double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2); // probably constant
      // (c) prior for alpha
      prod_k += gammapdfln(alpha(k), lambda_1, lambda_2);
    }
    int v = 0, prod_v = 0.0; // calc_loglik_prod_v();
    for(size_t v_sizet = 0; v_sizet < num_vocab; v_sizet++){
      v = static_cast<int>(v_sizet);
      for (int k = 0; k < k_seeded; k++){ // <--------- NOTE: using k_seeded
        // (a), seed topic part, second part numerator
        if (phi_s[k].find(v) == phi_s[k].end()){ // not a seed word
          prod_v += lgamma(beta_s);
        } else {
          prod_v += lgamma(beta_s + (double)n_x1_kv.coeffRef(k, v) );
        }
        // (a) regular part numerator
        prod_v += lgamma(beta + (double)n_x0_kv.coeffRef(k, v) );
      }
    }
    double prod_d = 0.0; // calc_loglik_prod_d();
    for (int d = 0; d < num_doc; d++){ // (c)
      prod_d += lgamma(alpha.sum());
      for(int k = 0; k < k_seeded; k++){ // <--------- NOTE: using k_seeded
        prod_d += lgamma(n_dk.coeffRef(d,k) + alpha(k)); // second numerator
        prod_d -= lgamma(alpha(k)); // first denominator
      }
      prod_d -= lgamma(n_dk.row(d).sum() + alpha.sum()); // second denominator
    }
    // (b) first part
    double others = (double)num_topics *
      ( lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) - lgamma(gamma_2)); // constant?
    double loglik = prod_k + prod_v + prod_d + others;
    double perplexity = exp(-loglik / (double)total_words);

    Rcerr << "log likelihood: " << loglik <<
             " (perplexity: " << perplexity << ")" << std::endl;

    checkUserInterrupt();
  }
  return model;
}

