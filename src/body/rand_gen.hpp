#pragma once
# include <iostream>
# include <random>
#include <vector>
#include <chrono>
using namespace std;

namespace randgen{
  // Default seed 
  int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937 mt(seed);
  mt19937 rand_gen = mt;

  void set_seed(int use_seed){
		static int set = 1;

		if(set){
			mt19937 mt(use_seed);
			seed = use_seed;
			rand_gen = mt;

			set = 0; // no more change
		}
  }
  int bernoulli(double &p){
    uniform_real_distribution<double> rand(0.0, 1.0);
    double r = rand(rand_gen);
    if(r <= p){
      return 1;
    }
    return 0;
  }
  double gamma(double a, double b){
    gamma_distribution<double> rand(a, 1.0 / b);
    return rand(rand_gen);
  }
  double beta(double a, double b){
    double ga = gamma(a, 1.0);
    double gb = gamma(b, 1.0);
    return ga / (ga + gb);
  }
  double uniform(double min = 0.0, double max = 1.0){
    uniform_real_distribution<double> rand(min, max);
    return rand(rand_gen);
  }
	int uniformint(int &size){
		uniform_int_distribution<int> rand(0, size);
		return rand(rand_gen);
	}
	vector<int> randomized_id_vec(int &num){
		vector<int> v(num);
		iota(v.begin(), v.end(), 0); // fill the vector
		shuffle(v.begin(), v.end(), rand_gen);
		return v;
	}
}

