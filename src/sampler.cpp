#include "sampler.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

namespace sampler{
	// inline int rand_wrapper(const int n) { return floor(unif_rand() * n); }

	double slice_uniform(double& lower, double& upper){
		return lower + (upper - lower) * unif_rand();
	}
}
