#ifndef __sampler__INCLUDED__
#define __sampler__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>

namespace sampler{
	inline int rand_wrapper(const int n) { return floor(unif_rand() * n); }
	double slice_uniform(double& lower, double& upper);
}

#endif
