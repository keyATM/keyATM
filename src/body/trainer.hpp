#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>

#include <boost/math/special_functions/digamma.hpp> // for digamma
#include <boost/math/special_functions/gamma.hpp> // for gamma functions
#include "rand_gen.hpp"

using namespace boost::math;

using namespace std;
using namespace Rcpp;
using namespace Eigen;



class Trainer{
public:
	FilesList *files_c;
	Vocab *vocab_c;

	// Data
	vector<vector<size_t>> W; // W_{d,i} Word IDs for each document
	unordered_map<string, int> doc_filename_to_id;
	unordered_map<int, string> doc_id_to_filename;
	vector<int> each_doc_len; // length of each document
	int total_words=0; // total number of words in data
	vector<unordered_map<size_t, double>> phi_s; // seed
	unordered_set<size_t> seed_words{}; // a vector of seed words
	vector<int> seed_num; // number of seeds for each topic
	int num_doc = 0; // D
	int num_vocab = 0; // V
	int num_seed_vocab = 0; // unique number of seeds
	int num_topics; // topic id starts from 0

	// Parameters
	vector<vector<int>> Z;
	vector<vector<int>> X;
	VectorXd alpha;
	double beta = 0.01; // currently, it's fixed
	double beta_s = 0.1;
	SparseMatrix<int, RowMajor> n_x0_kv; // K \times V sparse matrix, # of word v appear in regular topic k
	SparseMatrix<int, RowMajor> n_x1_kv; // K \times V sparse matrix, # of word v appear in seed topic k
	VectorXi n_x0_k; // K \tiems 1 vector, total number of words appear in regular topic k
	VectorXi n_x1_k; // K \tiems 1 vector, total number of words appear in seed topic k
	MatrixXd n_dk; // D \times K matrix, how many times topic k appear in a document d
	MatrixXd theta_dk; // theta, recoved for sampling alpha
	double gamma_1 = 1.0;
	double gamma_2 = 1.0;
	double lambda_1 = 1.0;
	double lambda_2 = 2.0;
	double seed_prob = 0.7; // for initialization

	// Tracking
	vector<int> iter_vec;
	vector<double> loglik_vec;
	vector<double> perplexity_vec;

	// Use During iteration
	int focus_docid;
	int focus_x;
	int focus_z;
	int focus_wid; // focusing word id
	size_t focus_wid_ui;
	int new_x;
	int new_z;
	double current_loglik=0.0;

	Trainer(string &datafolder, string &seed_path, int &rand_seed, int &iter_num);
	~Trainer();

	// Read Data
	void read_data();
	void read_doc(string &file_path);
	void add_sentence_to_doc(vector<string> &words, int doc_id);

	// Read Seeds
	void read_seed(string &seed_path);

	// Summary Output
	void read_data_summary();

	// Initialize
	void initialize();
	void initialize_alpha();
	void initialize_Z();
	void initialize_X();
	void initialize_counters();
	void count_n_x_kv(int &x, int &z, int &word_id);
	void count_n_x_k(int &x, int &z);
	void count_n_dk(int &doc_id, int &z);

	// Iteration
	void iteration(int &iter);
	void store_p_k(int &iter);
	void store_n_x0_k(int &iter);
	void store_n_x1_k(int &iter);
	void store_alpha(int &iter);
	// Remove data
	void remove_data();
	void remove_n_x_kv();
	void remove_n_x_k();
	void remove_n_dk();

	// Add Data
	void add_data();
	void add_n_x_kv();
	void add_n_x_k();
	void add_n_dk();

	// Update Z
	void update_Z(int &doc_id, int &w_position);
	int update_Z_x0();
	VectorXd calc_z_x0();
	int update_Z_x1();
	VectorXd calc_z_x1();

	// Update X
	void update_X(int &doc_id, int &w_position);
	double newprob_x0();
	double newprob_x1();

	// Update alpha
	void recover_theta();
	void update_alpha();
	double alpha_loglik(VectorXd &input_alpha);
	void slice_sampling_alpha(double min_v, double max_v, int max_shrink_time);


	// Calc_loglik and perplexity
	void tracking(int &iter);
	double calc_loglik();
	double calc_loglik_prod_k();
	double calc_loglik_prod_v();
	double calc_loglik_prod_d();
	double calc_loglik_others();
	double calc_perplexity(double &loglik);

	// Make output
	List get_output(int &show_words_num, bool &full_output);
	void make_Z_rlist(List &Z_rlist);
	void make_X_rlist(List &X_rlist);
	void make_TopWords(List &TopWords, int &show_words_num);
	void make_wordscount(vector< vector< pair<string,int> > > &wordscount);
	void sort_wordscount(vector< vector< pair<string,int> > > &wordscount);
	void make_TopWordsRList(vector< vector< pair<string,int> > > &wordscount,
													List &TopWords, int &show_words_num);
	void make_Words(List &RawWords, List &WordIDs);
	void make_iter_info(List &IterInfo);
};


