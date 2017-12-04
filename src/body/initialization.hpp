#pragma once
using namespace std;
#include "trainer.hpp"

void Trainer::initialize(){
	// alpha
	initialize_alpha();

	// Z
	initialize_Z();

	// X
	initialize_X();

	// n_x0_kv, n_x1_kv, n_x0_k, n_x1_k and n_dk
	initialize_counters();
}

void Trainer::initialize_alpha(){
	double alpha_k = 50.0 / num_topics;
	alpha = VectorXd::Constant(num_topics, alpha_k);
}

void Trainer::initialize_Z(){
	// Create a topic vector
	vector<int> topic_vec;
	for(int i=0; i < num_topics; i++ ){
		topic_vec.push_back(i);
	}

	// Randomly pick a topic
	int topic_id;
	for(size_t doc_id=0; doc_id<W.size(); ++doc_id){
		Z.push_back(vector<int>()); // add a vector for each document

		for(size_t word_position=0; word_position<W[doc_id].size(); ++word_position){
			topic_id = random_vec(topic_vec);
			Z[doc_id].push_back( topic_id );
		}
	}
}

void Trainer::initialize_X(){
	int x;
	for(size_t doc_id=0; doc_id<W.size(); ++doc_id){
		X.push_back(vector<int>()); // add a vector for each document

		for(size_t word_position=0; word_position<W[doc_id].size(); ++word_position){
			size_t word_id = W[doc_id][word_position];

			if(seed_words.find(word_id) == seed_words.end()){
				// Words are not in seed topics
				X[doc_id].push_back( 0 ); // words should be in regular topics
				continue;
			}else{
				x = regular_or_seed(seed_prob);
				X[doc_id].push_back( x );
			}
		}
	}
}

void Trainer::initialize_counters(){
	n_x0_kv = SparseMatrix<int, RowMajor> (num_topics, num_vocab);
	n_x1_kv = SparseMatrix<int, RowMajor> (num_topics, num_vocab);
	n_x0_k = VectorXi::Zero(num_topics);
	n_x1_k = VectorXi::Zero(num_topics);
	n_dk = MatrixXd::Zero(num_doc, num_topics);
	theta_dk = MatrixXd::Zero(num_doc, num_topics);

	int x;
	int z;
	int word_id;

	for(int doc_id=0; doc_id<num_doc; ++doc_id){
		for(int w_position=0; w_position<each_doc_len[doc_id]; ++w_position){
			x = X[doc_id][w_position];
			z = Z[doc_id][w_position];
			word_id = static_cast<int>(W[doc_id][w_position]);

			count_n_x_kv(x, z, word_id);
			count_n_x_k(x, z);
			count_n_dk(doc_id, z);
		}
	}

}

void Trainer::count_n_x_kv(int &x, int &z, int &word_id){
	if(x==0){
		n_x0_kv.coeffRef(z, word_id) += 1;
	}else if (x==1){
		n_x1_kv.coeffRef(z, word_id) += 1;
	}else{
		Rcout << "Something Goes Wrong" << endl;
	}
}

void Trainer::count_n_x_k(int &x, int &z){
	if(x==0){
		n_x0_k(z) += 1;
	}else if (x==1){
		n_x1_k(z) += 1;
	}else{
		Rcout << "Something Goes Wrong" << endl;
	}
}

void Trainer::count_n_dk(int &doc_id, int &z){
	n_dk.coeffRef(doc_id, z) += 1.0;
}


