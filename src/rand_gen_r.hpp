#pragma once
#include <iostream>
#include <random>
#include <vector>
#include <chrono>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

namespace randgen {
  // Default seed
  int seed = chrono::system_clock::now().time_since_epoch().count();

  void set_seed(int use_seed){
    // do nothing because we're set from R
  }

  int bernoulli(double & p){
    return R::unif_rand() <= p;
  }

  double gamma(double a, double b){
    return R::rgamma(a, 1.0 / b);
  }

  double beta(double a, double b){
    return R::rbeta(a, b);
  }

  double uniform(double min = 0.0, double max = 1.0){
    return R::runif(min, max);
  }

  // between 0 and size
  double uniformint(int & size){
    return floor(R::runif(0, size));
  }

  vector<int> randomized_id_vec(int & num){
    IntegerVector v = sample(num, num, false);
    vector<int> cpp_v = as<vector<int> >(v);
    return cpp_v;
  }
}
