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
//'
//' @keywords internal
// [[Rcpp::export]]
List read_dfm_cpp(Eigen::SparseMatrix<int> dfm, 
                  List W_raw, CharacterVector vocab, bool show_progress_bar)
{
  dfm = dfm.transpose();  // SparseMatrix is colmajor
  int doc_num = dfm.cols();
  string word_id;
  int count;

  Progress progress_bar(doc_num, show_progress_bar);

  for (int doc_id = 0; doc_id < doc_num; ++doc_id) {
    vector<string> doc_words;
    for (SparseMatrix<int>::InnerIterator it(dfm, doc_id); it; ++it) {
      word_id = vocab[it.index()];
      count = it.value();

      for (int i = 0; i < count; ++i) {
        doc_words.push_back(word_id);
      }
    }

    CharacterVector R_doc_words = Rcpp::wrap(doc_words);
    W_raw.push_back(R_doc_words);

    // Progress bar
    progress_bar.increment();
  }

  return W_raw;
}

