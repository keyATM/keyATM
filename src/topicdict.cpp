#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <string>

using namespace std;
using namespace Rcpp;
using namespace Eigen;

// Use c++11 and link functions to R
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppEigen")]]
// [[Rcpp::depends("BH")]]

#include "topicdict.hpp"

// [[Rcpp::export]]
List RunSeededLDA(std::string datafolder,
                  std::string seed_path, //changed to void from List
		              int iter_num=10,
		              int show_words_num=10,
		              bool full_output=true,
		              int seed=0){
	Rcout << "Read Data in " << datafolder << endl;
	Rcout << "Use seeds in " << seed_path << endl;

	Trainer trainer(datafolder, seed_path, seed, iter_num);

	// Initialize
	trainer.initialize();

	// Iteration
	Rcout << "Iteration:" << endl;
	//// Initial Loglik
	int first = 0;
	trainer.tracking(first);
	Rcout << endl;

	auto start = std::chrono::system_clock::now();
	for(int iter=0; iter<iter_num; ++iter){
		trainer.iteration(iter);

		if(iter % 10 == 0 && iter != 0){
			trainer.tracking(iter);
			auto dur = std::chrono::system_clock::now() - start;
			auto msec = std::chrono::duration_cast<std::chrono::seconds>(dur).count();
			Rcout << ", Time: " << msec << " sec" << endl;
			start = std::chrono::system_clock::now();
		}
	}

	// Make outputs
	List return_list = trainer.get_output(show_words_num, full_output);
	return return_list;
}
