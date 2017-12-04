#pragma once
using namespace std;
#include "trainer.hpp"

void Trainer::tracking(int &iter){
	iter_vec.push_back(iter);

	// double loglik = calc_loglik();
	// loglik_vec.push_back(loglik);

	////////////////////////////////////////
	// Here is an experimental code
	double loglik = 0.0;

	// using the polya distribution
	for(int k=0; k<num_topics; k++){

		for(int v=0; v<num_vocab; v++){
			// word
			loglik += lgamma(beta + (double)n_x0_kv.coeffRef(k, v) ) - lgamma(beta);
			loglik += lgamma(beta_s + (double)n_x1_kv.coeffRef(k, v) ) - lgamma(beta_s);
		}

		// word normalization
		loglik += lgamma( beta * (double)num_vocab ) - lgamma(beta * (double)num_vocab + (double)n_x0_kv.row(k).sum() );
		loglik += lgamma( beta_s * (double)num_vocab ) - lgamma(beta_s * (double)num_vocab + (double)n_x1_kv.row(k).sum() );

		// x
		loglik += lgamma( (double)n_x0_k(k) + gamma_2 ) - lgamma((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2)
				+ lgamma( (double)n_x1_k(k) + gamma_1 ) ;
		// x normalization
		loglik += lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) - lgamma(gamma_2);

	}

	// z
	for(int d=0; d<num_doc; d++){
		loglik += lgamma( alpha.sum() ) - lgamma( n_dk.row(d).sum() + alpha.sum() );

		for(int k=0; k<num_topics; k++){
			loglik += lgamma( n_dk.coeffRef(d,k) + alpha(k) ) - lgamma( alpha(k) );
		}
	}


	// using the "collapsed version": see Sato book page 127
	// double loglik_collapsed = 0.0;
	// double temp = 0.0;
	// int d = 0;
	// for(auto doc : W){
	// 	for(size_t word_id : doc){
	// 		int v = static_cast<int>(word_id);
	// 		temp = 0.0;

	// 		for(int k=0; k<num_topics; k++){
	// 			if(phi_s[k].find(v) == phi_s[k].end()){
	// 				// If a word is not seed
	// 				temp += ( ( n_dk.coeffRef(d,k) + alpha(k) ) / ( n_dk.row(k).sum() + alpha.sum() )     )  *
	// 					( beta + (double)n_x0_kv.coeffRef(k, v) ) / ( beta * (double)num_vocab + (double)n_x0_kv.row(k).sum() ) ;
	// 			}else{
	// 				// If a word is seed
	// 				temp += ( ( n_dk.coeffRef(d,k) + alpha(k) ) / ( n_dk.row(k).sum() + alpha.sum() )     )  *
	// 					( ( beta_s + (double)n_x1_kv.coeffRef(k, v) ) / ( beta_s * (double)num_vocab + (double)n_x1_kv.row(k).sum() ) *
	// 					( (double)n_x1_k(k) + gamma_1 ) / ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2) +
	// 					( beta + (double)n_x0_kv.coeffRef(k, v) ) / ( beta * (double)num_vocab + (double)n_x0_kv.row(k).sum() ) *
	// 					( (double)n_x0_k(k) + gamma_2 ) / ((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2)  ) ;

	// 			}

	// 		}

	// 		loglik_collapsed += log(temp);

	// 	}

	// 	++d;
	// }


	loglik_vec.push_back(loglik);
	////////////////////////////////////////

	double perplexity = calc_perplexity(loglik);
	perplexity_vec.push_back(perplexity);

	Rcout << "  Iteration Num: " << iter << ", Loglik/Perplexity: " << loglik << "/" << perplexity;
	// double perplexity_collapsed = calc_perplexity(loglik_collapsed);
	// cout << " collapsed: " << perplexity_collapsed;
}

double Trainer::calc_loglik(){
	// Check posterior distribution in derivation notes
	// There are three parts (a), (b) and (c)

	double prod_k = calc_loglik_prod_k();
	double prod_v = calc_loglik_prod_v();
	double prod_d = calc_loglik_prod_d();
	double others = calc_loglik_others();

	current_loglik = prod_k + prod_v + prod_d + others;
	Rcout << prod_k << " / " << prod_v << " / " << prod_d << " / " << others << "  ";
	return current_loglik;
}
double Trainer::calc_loglik_prod_k(){
	double sum = 0.0;

	for(int k=0; k<num_topics; ++k){
		// (a)
			// Seed Topic Part
			sum += lgamma((double)num_vocab * beta_s); // first term numerator
			sum -= lgamma(beta_s) * (double)num_vocab; // first term denominator
			sum -= lgamma( (double)num_vocab * beta_s + (double)n_x1_kv.row(k).sum()  ); // second term denominator

			// Regular Topic Part
			sum += lgamma((double)num_vocab * beta); //probably constant
			sum -= lgamma(beta) * (double)num_vocab; // probably constant
			sum -= lgamma( (double)num_vocab * beta + (double)n_x0_kv.row(k).sum()  ); // last part denominator

		// (b) second part
		sum += lgamma((double)n_x1_k(k) + gamma_1) + lgamma((double)n_x0_k(k) + gamma_2);
		sum -=  lgamma( (double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2); // probably constant

		// (c) prior for alpha
		sum += gammapdfln(alpha(k), lambda_1, lambda_2);
	}
	return sum;
}
double Trainer::calc_loglik_prod_v(){
	double sum = 0.0;
	int v = 0;

	for(size_t v_sizet=0; v_sizet<num_vocab; ++v_sizet){
		v = static_cast<int>(v_sizet);

		for(int k=0; k<num_topics; ++k){
			// (a), seed topic part, second part numerator
			if(phi_s[k].find(v) == phi_s[k].end()){
				// This is not a seed word
				sum += lgamma(beta_s);
			}else{
				// This is a seed word
				sum += lgamma(beta_s + (double)n_x1_kv.coeffRef(k, v) );
			}

			// (a) regular part numerator
			//cout <<  (double)n_x0_kv.coeffRef(k, v) << "/" << (double)n_x1_kv.coeffRef(k, v) << " " ;
			sum += lgamma(beta + (double)n_x0_kv.coeffRef(k, v) );
		}

	}
	return sum;

}
double Trainer::calc_loglik_prod_d(){
	double sum = 0.0;

	// (c)
	for(int d=0; d<num_doc; ++d){
		sum += lgamma(alpha.sum());

		for(int k=0; k<num_topics; ++k){
			sum += lgamma(n_dk.coeffRef(d,k) + alpha(k));	 // second numerator
			sum -= lgamma(alpha(k)); // first denominator

		}

		sum -= lgamma(n_dk.row(d).sum() + alpha.sum()); // second denominator

	}
	return sum;
}
double Trainer::calc_loglik_others(){
	double sum = 0.0;

	// (b) first part
	sum += (double)num_topics * ( lgamma(gamma_1 + gamma_2) - lgamma(gamma_1) -  lgamma(gamma_2)); //probably constant

	return sum;
}
double Trainer::calc_perplexity(double &loglik){
	double perplexity = 0.0;

	perplexity = exp(-loglik / (double)total_words);

	return perplexity;
}


