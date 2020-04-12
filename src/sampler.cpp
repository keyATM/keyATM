#include "sampler.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;

namespace sampler{
	double u;
	double temp;
	int index;

	double slice_uniform(double& lower, double& upper)
  {
		return lower + (upper - lower) * R::unif_rand();
	}


	std::vector<int> shuffled_indexes(int m)
  {
    // Returns a vector of shuffled indexes for sampling
		std::vector<int> v(m);
		std::iota(v.begin(), v.end(), 0);
		std::random_shuffle(v.begin(), v.end(), sampler::rand_wrapper);
		return v;
	}


	int rcat(Eigen::VectorXd &prob, int &size)
  { 
		u = R::runif(0, 1);
		temp = 0.0;
		index = 0;
		for (int ii = 0; ii < size; ii++) {
			temp += prob(ii);
			if (u < temp) {
				index = ii;
				break;
			}
		}
		return index;
	}


	int rcat_without_normalize(Eigen::VectorXd &prob, double &total, int &size)
  { 
    // Draw from a categorial distribution
    // This function does not requiare a normalized probability vector.
		u = R::runif(0, 1) * total;
		temp = 0.0;
		index = 0;
		for (int ii = 0; ii < size; ii++) {
			temp += prob(ii);
			if (u < temp) {
				index = ii;
				break;
			}
		}
		return index;
	}
}
