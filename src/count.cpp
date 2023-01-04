#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

// Use c++11 and link functions to R
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
bool word_in_doc(StringVector doc, std::string word)
{
  int size = doc.length();
  for (int i = 0; i < size; ++i) {
    if (doc[i] == word) {
      return true;
    }
  }
  return false;
}
