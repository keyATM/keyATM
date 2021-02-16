#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <progress.hpp>
#include <progress_bar.hpp>


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Eigen;
using namespace Rcpp;
using namespace std;


//' Read files from the quanteda dfm (this is the same as dgCMatrix) 
//'
//' @param dfm a dfm input (sparse Matrix)
//' @param W_raw an object to return
//' @param vocab a vector of vocabulary
//' @param show_progress_bar show a progress bar 
//' @param split split proportion 
//'
//' @keywords internal
// [[Rcpp::export]]
List read_dfm_cpp(Eigen::SparseMatrix<int> dfm, 
                  List W_read, CharacterVector vocab, bool show_progress_bar,
                  double split)
{
  dfm = dfm.transpose();  // SparseMatrix is colmajor
  int doc_num = dfm.cols();
  string word_id;
  int count;
  double u;

  List W_raw = W_read["W_raw"];
  List W_split = W_read["W_split"];

  Progress progress_bar(doc_num, show_progress_bar);

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

    CharacterVector R_doc_words = Rcpp::wrap(doc_words);
    W_raw.push_back(R_doc_words);

    if (split != 0) {
      CharacterVector R_doc_words_split = Rcpp::wrap(doc_words_split); 
      W_split.push_back(R_doc_words_split);
    }

    // Progress bar
    progress_bar.increment();
  }

  // make an output
  W_read["W_raw"] = W_raw;
  W_read["W_split"] = W_split;

  return W_read;
}

