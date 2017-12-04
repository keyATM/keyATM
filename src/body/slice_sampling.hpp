#pragma once
#include <functional>
#include <math.h>
#include <algorithm>
using namespace std;
using namespace Eigen;

#include "trainer.hpp"

double expandp(double p){
    // p --> (x > 0)
    return p / (1.0 - p);
}

double shrinkp(double x){
    // (x > 0) --> p
    return x / (1.0 + x);
}

void Trainer::slice_sampling_alpha(double min_v = 1e-9, double max_v=20.0, int max_shrink_time=3000)
{
	int k;
	double start;
	double end;
	double previous_p;
	double new_p;
	double newlikelihood;
	double slice_;
	VectorXd keep_current_param = alpha;

  vector<int> topic_ids(num_topics);
  iota(topic_ids.begin(), topic_ids.end(), 0); // 0: starting number, fill the vector
	std::shuffle(topic_ids.begin(), topic_ids.end(), randgen::rand_gen);


	for(int i=0; i<num_topics; i++){
		k = topic_ids[i];

		start = shrinkp(min_v);
		end = 1.0;
		// end = shrinkp(max_v);

		previous_p = shrinkp(alpha(k));

		slice_ = alpha_loglik(alpha) - 2.0 * log(1.0 - previous_p) + log(randgen::uniform());

		for(int shrink_time=0; shrink_time<max_shrink_time; shrink_time++){ // for each dimension
			new_p = randgen::uniform(start, end);
			alpha(k) = expandp(new_p);

			newlikelihood = alpha_loglik(alpha) - 2.0 * log(1.0 - new_p);

			if (slice_ < newlikelihood){
				break;
			}else if(previous_p < new_p){
				end = new_p;
			}else if(new_p < previous_p){
				start = new_p;
			}else{
				alpha(k) = keep_current_param(k);
			}

		}

	}

}
