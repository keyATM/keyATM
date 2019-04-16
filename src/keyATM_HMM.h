#ifndef __keyATM_HMM__INCLUDED__
#define __keyATM_HMM__INCLUDED__

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include "sampler.h"
#include "keyATM.h"

using namespace Eigen;
using namespace Rcpp;
using namespace std;


class keyATMhmm : public keyATMbase
{
	public:

		// Constructor
		keyATMhmm(List model_, const int iter_, const int output_per_);

};


#endif
