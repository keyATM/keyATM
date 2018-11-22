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
}
