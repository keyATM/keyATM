#include "sampler.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

namespace sampler{
	double slice_uniform(double& lower, double& upper){
		return lower + (upper - lower) * unif_rand();
	}


	std::vector<int> shuffled_indexes(int m){
		std::vector<int> v(m);
		std::iota(v.begin(), v.end(), 0);
		std::random_shuffle(v.begin(), v.end(), sampler::rand_wrapper);
		return v;
	}

	int rcat(Eigen::VectorXd &prob){ // was called 'multi1'
		double u = R::runif(0, 1);
		double temp = 0.0;
		int index = 0;
		for (int ii = 0; ii < prob.size(); ii++){
			temp += prob(ii);
			if (u < temp){
				index = ii;
				break;
			}
		}
		return index;
	}

	int rcat_without_normalize(Eigen::VectorXd &prob, double &total){ // was called 'multi1'
		double u = R::runif(0, 1) * total;
		double temp = 0.0;
		int index = 0;
		for (int ii = 0; ii < prob.size(); ii++){
			temp += prob(ii);
			if (u < temp){
				index = ii;
				break;
			}
		}
		return index;
	}
}
