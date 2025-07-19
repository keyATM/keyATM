#include <Rcpp.h>

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>

#include <algorithm>
#include <iostream>
#include <string>

// Progress bar
#include <cli/progress.h>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace Rcpp;
using namespace std;

//' Read files from the quanteda dfm (this is the same as dgCMatrix)
//'
//' @param dfm a dfm input (sparse Matrix)
//' @param W_read an object to return
//' @param vocab a vector of vocabulary
//' @param split split proportion
//'
//' @keywords internal
// [[Rcpp::export]]
List read_dfm_cpp(Eigen::SparseMatrix<int> dfm, List W_read, StringVector vocab,
                  double split) {
  dfm = dfm.transpose(); // SparseMatrix is colmajor
  int doc_num = dfm.cols();
  string word_id;
  int count;
  double u;

  List W_raw = W_read["W_raw"];
  List W_split = W_read["W_split"];

  SEXP progress_bar = PROTECT(cli_progress_bar(doc_num, NULL));
  cli_progress_set_name(progress_bar, "Loading documents");

  for (int doc_id = 0; doc_id < doc_num; ++doc_id) {
    vector<string> doc_words;
    vector<string> doc_words_split;

    for (SparseMatrix<int>::InnerIterator it(dfm, doc_id); it; ++it) {
      word_id = vocab[it.index()];
      count = it.value();

      if (split != 0) {
        // Split the dfm

        for (int i = 0; i < count; ++i) {
          u = R::unif_rand();

          if (u < split) {
            // Smaller subset
            doc_words_split.push_back(word_id);
          } else {
            doc_words.push_back(word_id);
          }
        }
      } else {
        // Do not split the dfm
        for (int i = 0; i < count; ++i) {
          doc_words.push_back(word_id);
        }
      }
    }

    StringVector R_doc_words = Rcpp::wrap(doc_words);
    W_raw.push_back(R_doc_words);

    if (split != 0) {
      StringVector R_doc_words_split = Rcpp::wrap(doc_words_split);
      W_split.push_back(R_doc_words_split);
    }

    if (CLI_SHOULD_TICK) {
      cli_progress_set(progress_bar, doc_id);
    }

    // Check keybord interruption to cancel the iteration
    checkUserInterrupt();
  }

  cli_progress_done(progress_bar);
  UNPROTECT(1);

  // make an output
  W_read["W_raw"] = W_raw;
  W_read["W_split"] = W_split;

  return W_read;
}
