#pragma once
#include <functional>
using namespace std;
#include "trainer.hpp"
#include "slice_sampling.hpp"
#include <math.h>

void Trainer::iteration(int &iter){

	int doc_id;
	int w_position;
	vector<int> randomized_docid = randomized_id_vec(num_doc);
	vector<int> randomized_wordpos;

	for(int doc_index=0; doc_index<num_doc; ++doc_index){
		doc_id = randomized_docid[doc_index];
		randomized_wordpos = randomized_id_vec(each_doc_len[doc_id]);

		for(int w_index=0; w_index<each_doc_len[doc_id]; ++w_index){
			w_position = randomized_wordpos[w_index];
			update_Z(doc_id, w_position);
			update_X(doc_id, w_position);
		}
	}

	// Hyper parameter
	update_alpha();

	// Debug
	store_p_k(iter);
	store_n_x0_k(iter);
	store_n_x1_k(iter);
	store_alpha(iter);
}

void Trainer::store_p_k(int &iter){
	VectorXd p_k = VectorXd::Zero(num_topics);
	for(int k=0; k<num_topics; k++){
		p_k(k) = (double)n_x1_k(k) / ((double)n_x0_k(k) + (double)n_x1_k(k));
	}

	iter_log::store_p_k(iter, p_k);
}

void Trainer::store_alpha(int &iter){
	iter_log::store_alpha(iter, alpha);
}

void Trainer::store_n_x0_k(int &iter){
	iter_log::store_n_x0_k(iter, n_x0_k);
}

void Trainer::store_n_x1_k(int &iter){
	iter_log::store_n_x1_k(iter, n_x1_k);
}

// Remove data
void Trainer::remove_data(){
	remove_n_x_kv();
	remove_n_x_k();
	remove_n_dk();
}
void Trainer::remove_n_x_kv(){
	if(focus_x==0){
		n_x0_kv.coeffRef(focus_z, focus_wid) -= 1;
	}else if(focus_x==1){
		n_x1_kv.coeffRef(focus_z, focus_wid) -= 1;
	}
}
void Trainer::remove_n_x_k(){
	if(focus_x==0){
		n_x0_k(focus_z) -= 1;
	}else if(focus_x==1){
		n_x1_k(focus_z) -= 1;
	}
}
void Trainer::remove_n_dk(){
	n_dk.coeffRef(focus_docid, focus_z) -= 1;
}

// Add Data
void Trainer::add_data(){
	add_n_x_kv();
	add_n_x_k();
	add_n_dk();
}
void Trainer::add_n_x_kv(){
	if(focus_x==0){
		n_x0_kv.coeffRef(focus_z, focus_wid) += 1;
	}else if(focus_x==1){
		n_x1_kv.coeffRef(focus_z, focus_wid) += 1;
	}
}
void Trainer::add_n_x_k(){
	if(focus_x==0){
		n_x0_k(focus_z) += 1;
	}else if(focus_x==1){
		n_x1_k(focus_z) += 1;
	}
}
void Trainer::add_n_dk(){
	n_dk.coeffRef(focus_docid, focus_z) += 1;
}

// Update Z
void Trainer::update_Z(int &doc_id, int &w_position){
	focus_docid = doc_id;
	focus_x = X[doc_id][w_position];
	focus_z = Z[doc_id][w_position];
	focus_wid = static_cast<int>(W[doc_id][w_position]);
	focus_wid_ui = W[doc_id][w_position];
	remove_data();

	if(focus_x==0){
		new_z = update_Z_x0();
	}else if(focus_x==1){
		new_z = update_Z_x1();
	}else{
		Rcout << "Something goes wrong" << endl;
	}

	// Update
	Z[doc_id][w_position] = new_z;
	focus_z = new_z;
	add_data();
}
int Trainer::update_Z_x0(){
	VectorXd z_prob_vec = calc_z_x0();
	return multi1(z_prob_vec);
}
VectorXd Trainer::calc_z_x0(){
	VectorXd z_prob_vec = VectorXd::Zero(num_topics);

	for(int k=0; k<num_topics; ++k){
		double numerator = log(beta + (double)n_x0_kv.coeffRef(k, focus_wid)) +
											log((double)n_x0_k(k) + gamma_2) +
											log((double)n_dk.coeffRef(focus_docid, k) + alpha(k));

		double denominator = log((double)num_vocab * beta + (double)n_x0_kv.row(k).sum() ) +
					log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);

		z_prob_vec(k) = numerator - denominator;
	}

	// Normalization
	double sum = logsumexp_Eigen(z_prob_vec);
	for(int k=0; k<num_topics; ++k)
		z_prob_vec(k) = exp( z_prob_vec(k) - sum );

	//cout << z_prob_vec << endl; // for debug

	return z_prob_vec;
}
int Trainer::update_Z_x1(){
	VectorXd z_prob_vec = calc_z_x1();
	return multi1(z_prob_vec);
}
VectorXd Trainer::calc_z_x1(){
	VectorXd z_prob_vec = VectorXd::Zero(num_topics);
	vector<int> make_zero_later;
	double numerator;
	double denominator;

	for(int k=0; k<num_topics; ++k){

		if(phi_s[k].find(focus_wid_ui) == phi_s[k].end()){
			z_prob_vec(k) = 1.0;
			make_zero_later.push_back(k);
			continue;
		}else{
			numerator = log(beta_s + (double)n_x1_kv.coeffRef(k, focus_wid)) +
							log( ((double)n_x1_k(k) + gamma_1) ) +
							log( ((double)n_dk.coeffRef(focus_docid, k) + alpha(k)) );
		}


		denominator = log((double)seed_num[k] * beta_s + (double)n_x1_kv.row(k).sum() ) +
			log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);

		z_prob_vec(k) = numerator - denominator;
	}

	// Normalization
	double sum = logsumexp_Eigen(z_prob_vec);
	for(int k=0; k<num_topics; ++k)
		z_prob_vec(k) = exp( z_prob_vec(k) - sum );

	//	Delete the commented out lines below later
	for(int k : make_zero_later)
		z_prob_vec(k) = 0.0;

	z_prob_vec = z_prob_vec / z_prob_vec.sum();

	//cout << z_prob_vec.transpose() << endl; // for debug

	return z_prob_vec;
}

// Update X
void Trainer::update_X(int &doc_id, int &w_position){
	focus_docid = doc_id;
	focus_x = X[doc_id][w_position];
	focus_z = Z[doc_id][w_position];
	focus_wid = static_cast<int>(W[doc_id][w_position]);
	focus_wid_ui = W[doc_id][w_position];
	remove_data();

	double x1_logprob = newprob_x1();

	if(x1_logprob == -1.0){
		// If the probability of x_di = 1 case is 0, it should be x=0 (regualr topic)
		new_x = 0;
	}else{
		double x0_logprob = newprob_x0();
		// Normalize
		double sum = 0.0;
		double prob_array[] = {x0_logprob, x1_logprob};
		for (int i = 0; i < 2; ++i)
			sum = logsumexp (sum, prob_array[i], (i == 0));

		double x0_prob = exp( x0_logprob - sum);
		double x1_prob = exp( x1_logprob - sum);
		new_x = Bern(x0_prob, x1_prob);

		//cout << "  " <<  x0_prob << " / " << x1_prob << " / " << new_x << endl; // for debug
	}

	// Update
	X[doc_id][w_position] = new_x;
	focus_x = new_x;
	add_data();

}
double Trainer::newprob_x0(){
	double prob;
	int k = focus_z;

	double numerator = log(beta + (double)n_x0_kv.coeffRef(k, focus_wid)) +
										log((double)n_x0_k(k) + gamma_2);


	double denominator = log((double)num_vocab * beta + (double)n_x0_kv.row(k).sum() ) +
				log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);

	//cout << numerator << " / " << denominator << endl;
	prob = numerator - denominator;


	return prob;
}
double Trainer::newprob_x1(){
	double prob;
	int k = focus_z;
	double numerator;
	double denominator;


	if(phi_s[k].find(focus_wid_ui) == phi_s[k].end()){
		return -1.0;
	}else{
		numerator = log(beta_s + (double)n_x1_kv.coeffRef(k, focus_wid)) +
						log( ((double)n_x1_k(k) + gamma_1) );
	}

	denominator = log((double)seed_num[k] * beta_s + (double)n_x1_kv.row(k).sum() ) +
		log((double)n_x1_k(k) + gamma_1 + (double)n_x0_k(k) + gamma_2);

	prob = numerator - denominator;

	return prob;
}

void Trainer::update_alpha(){
	slice_sampling_alpha();
}


double Trainer::alpha_loglik(VectorXd &input_alpha){
	// // Recover Theta
	// double loglik = 0.0;
	// double fixed = 0.0;
	// VectorXd theta_d;
	//
	//
	// for(int k=0; k<num_topics; k++){
	// 	// Add prior
	// 	loglik += gammapdfln(input_alpha(k), lambda_1, lambda_2);
	// }
	//
	//
	// for(int d=0; d<num_doc; d++){
	// 	fixed = lgamma( ( n_dk.row(d).array() + input_alpha.array()).sum() );
	// 	for(int k=0; k<num_topics; k++){
	// 		fixed -= lgamma( n_dk.coeffRef(d, k) + input_alpha(k) );
	// 	}
	//
	// 	theta_d = ( n_dk.row(d).array() + input_alpha.array() ) / ( (double)each_doc_len[d] + input_alpha.sum() );
	//
	// 	loglik += ( (n_dk.row(d).array() + input_alpha.array() -  1.0  ).array() * theta_dk.array() ).sum();
	// 	loglik += fixed;
	// }
	//
	// return loglik;

	// Collapsed
	double loglik = 0.0;
	double fixed_part = 0.0;
	VectorXd ndk_ak;

	fixed_part += lgamma( input_alpha.sum() ); // first term numerator

	for(int k=0; k<num_topics; k++){
		fixed_part -= lgamma( input_alpha(k) ); // first term denominator

		// Add prior
		loglik += gammapdfln(input_alpha(k), lambda_1, lambda_2);
	}

	for(int d=0; d<num_doc; d++){
		loglik += fixed_part;

		ndk_ak = n_dk.row(d).transpose() + input_alpha;

		// second term numerator
		for(int k=0; k<num_topics; k++){
			loglik += lgamma( ndk_ak(k) );
		}

		// second term denominator
		loglik -= lgamma( ndk_ak.sum() );
	}

	return loglik;

}
