#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Eigen;
using namespace Rcpp;


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
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

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double logsumexp_Eigen(VectorXd &vec){
  double sum = 0.0;
  int index;
  for(size_t i = 0; i < vec.size(); i++){
    index = static_cast<int>(i);
    sum = logsumexp(sum, vec[index], (index == 0));
  }
  return sum;
}


// [[Rcpp::depends(RcppEigen)]]
int rcat(VectorXd &prob){
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



// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List train(List model, int k_seeded, int k_free, double alpha_k,
                  int iter = 0){

  List W = model["W"], Z = model["Z"], X = model["X"];
  StringVector files = model["files"], vocab = model["vocab"];
  List seeds = model["seeds"]; // convert this to T&S's phi form

  std::vector< std::map<int, double> > phi_s(seeds.size()); // change to unordered map later
  std::vector<int> seed_num(seeds.size());
  for (int ii = 0; ii < seeds.size(); ii++){
    IntegerVector wd_ids = seeds[ii];
    seed_num[ii] = wd_ids.size();
    std::map<int, double> phi_sk;
    for (int jj = 0; jj < wd_ids.size(); jj++)
      phi_sk[wd_ids(jj)] = 1.0 / wd_ids.size();
    phi_s[ii] = phi_sk;
  }

  // alpha-related constants
  int num_topics = k_seeded + k_free;
  VectorXd alpha = VectorXd::Constant(num_topics, alpha_k);

  // document-related constants
  int num_vocab = vocab.size();
  int num_doc = files.size();

  // distributional constants
  double gamma_1 = 1.0, gamma_2 = 1.0;
  double lambda_1 = 1.0, lambda_2 = 2.0;
  double beta = 0.01, beta_s = 0.1;

  // storage for sufficient statistics
  SparseMatrix<int, RowMajor> n_x0_kv = SparseMatrix<int, RowMajor> (num_topics, num_vocab);
  SparseMatrix<int, RowMajor> n_x1_kv = SparseMatrix<int, RowMajor> (num_topics, num_vocab);
  MatrixXd n_dk = MatrixXd::Zero(num_doc, num_topics);
  MatrixXd theta_dk = MatrixXd::Zero(num_doc, num_topics); // TODO initialize?
  // and their margins
  VectorXi n_x0_k = VectorXi::Zero(num_topics);
  VectorXi n_x1_k = VectorXi::Zero(num_topics);

  // fill n_x0_kv, n_x0_k, n_dk
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

  // Note: no randomized update sequence in dev
  for (int ii = 0; ii < iter; ii++){
    for (int doc_id = 0; doc_id < num_doc; doc_id++){
      IntegerVector doc_x = X[doc_id], doc_z = Z[doc_id], doc_w = W[doc_id];
      for (int w_position = 0; w_position < doc_x.size(); w_position++){
        int x = doc_x[w_position], z = doc_z[w_position], w = doc_w[w_position];

        // remove data
        if (x == 0){
          n_x0_kv.coeffRef(z, w_position) -= 1;
          n_x0_k(z) -= 1;
        } else {
          n_x1_kv.coeffRef(z, w_position) -= 1;
          n_x1_k(z) -= 1;
        }
        n_dk.coeffRef(doc_id, z) -= 1;

        // draw z
        VectorXd z_prob_vec = VectorXd::Zero(num_topics);
        int new_z = -1; // debug
        if (x == 0){
          for (int k = 0; k < num_topics; k++){
            double numerator = log(beta + (double)n_x0_kv.coeffRef(k, w_position)) +
              log((double)n_x0_k(k) + gamma_2) +
              log((double)n_dk.coeffRef(doc_id, k) + alpha(k));
            double denominator = log((double)num_vocab * beta + (double)n_x0_kv.row(k).sum()) +
              log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);
            z_prob_vec(k) = numerator - denominator;
          }
          double sum = logsumexp_Eigen(z_prob_vec); // normalize
          for (int k = 0; k < num_topics; k++)
            z_prob_vec(k) = exp(z_prob_vec(k) - sum);
          new_z = rcat(z_prob_vec); // take a sample

        } else {
          std::vector<int> make_zero_later; // elements to zero out
          for (int k = 0; k < num_topics; k++){
            double numerator; // is this right? It will be zero, presumably
            if (phi_s[k].find(w_position) == phi_s[k].end()){ // w isn't a seed
              z_prob_vec(k) = 1.0;
              make_zero_later.push_back(k);
              continue;
            } else{ // w is a seed
              numerator = log(beta_s + (double)n_x1_kv.coeffRef(k, w_position)) +
                log( ((double)n_x1_k(k) + gamma_1) ) +
                log( ((double)n_dk.coeffRef(doc_id, k) + alpha(k)) );
            }
            double denominator = log((double)seed_num[k] * beta_s + (double)n_x1_kv.row(k).sum() ) +
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
        // swap in the new z
        doc_z[w_position] = new_z;
        z = new_z;

        // add back data counts
        if (x == 0){
          n_x0_kv.coeffRef(z, w_position) += 1;
          n_x0_k(z) += 1;
        } else {
          n_x1_kv.coeffRef(z, w_position) += 1;
          n_x1_k(z) += 1;
        }
        n_dk.coeffRef(doc_id, z) += 1;

        //update_X(doc_id, w_position);

      }
    }
    // update_alpha();
  }

  return model;
}


//vector<int> each_doc_len(num_doc);
//for(int doc_id = 0; doc_id < num_doc; doc_id++){
//  IntegerVector w = W[doc_id];
//  each_doc_len[doc_id] = w.size();
//}
