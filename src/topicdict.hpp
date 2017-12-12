#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/math/special_functions/digamma.hpp> // for digamma
#include <boost/math/special_functions/gamma.hpp> // for gamma functions
#include <dirent.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <math.h>
#include <numeric>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <random>
#include <chrono>

#include "rand_gen_r.hpp"
//#include "rand_gen.hpp" /* original c++ random number generators */

using namespace std;
using namespace Eigen;
using namespace Rcpp;
using namespace boost::math;

// [[Rcpp::plugins("cpp11")]]

/* LOGGING NAMESPACE */
namespace iter_log{
  NumericMatrix iter_alpha;
  NumericMatrix iter_n_x0_k;
  NumericMatrix iter_n_x1_k;
  NumericMatrix iter_p_k;

  template <class EigenVec>
  void store_values(NumericMatrix &rmat, int &iter_num, EigenVec &vec);

  void initialize(int &iter_num, int &num_topics){
    static int set = 1;

    if(set){
      NumericMatrix temp(iter_num, num_topics);
      iter_alpha = clone(temp);
      iter_n_x0_k = clone(temp);
      iter_n_x1_k = clone(temp);
      iter_p_k = clone(temp);

      set = 0; // initialize only once
    }
  }

  void store_alpha(int &iter_num, VectorXd &vec){
    store_values(iter_alpha, iter_num, vec);
  }

  void store_n_x0_k(int &iter_num, VectorXi &vec){
    store_values(iter_n_x0_k, iter_num, vec);
  }

  void store_n_x1_k(int &iter_num, VectorXi &vec){
    store_values(iter_n_x1_k, iter_num, vec);
  }

  void store_p_k(int &iter_num, VectorXd &vec){
    store_values(iter_p_k, iter_num, vec);
  }

  template <class EigenVec>
  void store_values(NumericMatrix &rmat, int &iter_num, EigenVec &vec){
    int num = vec.size();

    for(int i=0; i<num; i++){
      rmat(iter_num, i) = vec(i);
    }
  }
}


// Create Output
void output_clear(const string &file_name, int col_num=0){
  ofstream outfile(file_name, ios_base::out);

  if(col_num==0){
    outfile << "parameter,value"<<endl;
  }else{
    for(int c=0; c<col_num-1; c++){
      outfile << "Col" << c << ",";
    }
    outfile<< "Col" << col_num-1 << endl;
  }

  outfile.close();
}

void output_value(const string &file_name, double value, int iter){
  ofstream outfile(file_name, ios_base::app);

  outfile << iter << "," << value <<endl;
  outfile.close();
}

void output_EigenVec(const string &file_name, VectorXd &vec){ // for float value
  ofstream outfile(file_name, ios_base::app);

  IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ",", ",", "", "", "", "");
  outfile << vec.transpose().format(CommaInitFmt) <<endl;
  outfile.close();
}

void output_EigenVec(const string &file_name, VectorXi &vec){ // for int value
  ofstream outfile(file_name, ios_base::app);

  IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ",", ",", "", "", "", "");
  outfile << vec.transpose().format(CommaInitFmt) <<endl;
  outfile.close();
}

// Stats func
int random_vec(vector<int> &vec){
  int length = static_cast<int>( vec.size() - 1 );
  int index = randgen::uniformint(length);
  int random_element = vec[index];
  return random_element;
}

int regular_or_seed(double &prob){
  double value = randgen::uniform();

  if(value <= prob){
    return 1;
  }else{
    return 0;
  }
}

vector<int> randomized_id_vec(int &num){
  return randgen::randomized_id_vec(num);
}

int Bern(double &prob0, double &prob1){
  return randgen::bernoulli(prob1); /// !!!!
}

int multi1(VectorXd &prob){
  // Multi(x, 1), return category index
  double u = randgen::uniform();
  double temp = 0.0;
  int index = 0;
  for(int i=0; i<prob.size(); i++){
    temp += prob(i);
    if(u < temp){
      index = i;
      break;
    }
  }
  return index;
}

double gammapdfln(double x, double a, double b){
  return a * log(b) - lgamma(a) + (a-1.0) * log(x) - b * x;
}

// Strings
void split_string_by(const string &str, char delim, vector<string> &elems){
  elems.clear();
  string item;
  for(char ch: str){
    if (ch == delim){
      if (!item.empty()){
        elems.push_back(item);
      }
      item.clear();
    }
    else{
      item += ch;
    }
  }
  if (!item.empty()){
    elems.push_back(item);
  }
}

void split_word_by(const string &str, char delim, vector<string> &elems){
  elems.clear();
  string item;
  for(char ch: str){
    if (ch == delim){
      if (!item.empty()){
        elems.push_back(item);
      }
      item.clear();
    }
    else{
      item += ch;
    }
  }
  if (!item.empty()){
    elems.push_back(item);
  }
}


double logsumexp (double &x, double &y, bool flg){
  if (flg) return y; // init mode
  if (x == y) return x + 0.69314718055; // log(2)
  double vmin = std::min (x, y);
  double vmax = std::max (x, y);
  if (vmax > vmin + 50) {
    return vmax;
  } else {
    return vmax + std::log (std::exp (vmin - vmax) + 1.0);
  }
}

double logsumexp_Eigen(VectorXd &vec){
  double sum = 0.0;
  int index;
  for(size_t i=0; i<vec.size(); ++i){
    index = static_cast<int>(i);
    sum = logsumexp (sum, vec[index], (index == 0));
  }
  return sum;
}


// STRING COMPARISONS
bool sortbysec_descend(const pair<string,int> &a, const pair<string,int> &b){
  // Sort pair by the second element, descending
  return (a.second > b.second);
}

struct comp{
  comp(string const& s) : _s(s) { }

  bool operator () (pair<string,int> const& p){
    return (p.first == _s);
  }
  std::string _s;
};



class FilesList{
public:
	vector<string> files_list; // stores path to the files we want to use
	int file_num;

	FilesList(string &folder_path){
		files_list = getFileList(folder_path);
	}

	vector<string> getFileList(const string &folder_path) {
		// Get File list: http://xr0038.hatenadiary.jp/entry/2015/02/26/154608sR

	  DIR *dp;       // pointer to the directory
	  dirent* entry; // readdir() returns this entry point

	  dp = opendir(folder_path.c_str()); // don't forget to change string to const* char
		vector<string> path_list;
	  if (dp==NULL) {
		  Rcpp::stop("Folder does not exist");
	  }
	  do {
	    entry = readdir(dp);
	    if (entry != NULL){
				if(has_suffix(entry->d_name, ".txt"))
					path_list.push_back(folder_path + entry->d_name); // make full path
					++file_num;
			}
	  } while (entry != NULL);

		return path_list;
	}

	bool has_suffix(const string& s, const string& suffix){ // check suffix
			return (s.size() >= suffix.size()) && equal(suffix.rbegin(), suffix.rend(), s.rbegin());
	}

};


class Vocab{
public:
	unordered_map<size_t, string> wordid_to_word;
	unordered_map<string, size_t> word_to_wordid;
	int num_vocab = 0;

	Vocab(){}

	size_t get_wordid(string &word){
		if(word_to_wordid.find(word) == word_to_wordid.end()){
			// New word
			size_t word_id = word_to_wordid.size();
			word_to_wordid[word] = word_id;
			wordid_to_word[word_id] = word;

			++ num_vocab;
		}
		return word_to_wordid[word];
	}
	string get_word(size_t &word_id){
		return wordid_to_word[word_id];
	}
	int get_num_vocab(){
		return num_vocab;
	}
};

class Trainer{
public:
  FilesList *files_c;
  Vocab *vocab_c;
  //FilesList files_c;
  //Vocab vocab_c;



  // Data
  vector<vector<size_t> > W; // W_{d,i} Word IDs for each document
  unordered_map<string, int> doc_filename_to_id;
  unordered_map<int, string> doc_id_to_filename;
  vector<int> each_doc_len; // length of each document
  int total_words=0; // total number of words in data
  vector<unordered_map<size_t, double> > phi_s; // seed
  unordered_set<size_t> seed_words{}; // a vector of seed words
  vector<int> seed_num; // number of seeds for each topic
  int num_doc = 0; // D
  int num_vocab = 0; // V
  int num_seed_vocab = 0; // unique number of seeds
  int num_topics; // topic id starts from 0

  // Parameters
  vector<vector<int> > Z;
  vector<vector<int> > X;
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

// Constructor and Destructor
Trainer::Trainer(string &datafolder, string &seed_path, int &rand_seed, int &iter_num){ // constructor

  // Set seed
  if (rand_seed != 0){
    randgen::set_seed(rand_seed);
  }

  // Make a file list
  files_c = new FilesList(datafolder);

  // Read files
  vocab_c = new Vocab;
  read_data();

  // Read Seed
  read_seed(seed_path);

  // Summary Output
  read_data_summary();

  // Prepare logger
  iter_log::initialize(iter_num, num_topics);

}

Trainer::~Trainer(){
  delete files_c;
  delete vocab_c;
}



///////////////
// Read Data //
///////////////

void Trainer::read_data(){
	//cout << "Read following files:\n";
	for(auto file_path : files_c -> files_list){
		//cout << "  " << file_path << "\n";
		read_doc(file_path);
		++num_doc;
	}
	Rcout << endl;
	num_vocab = vocab_c -> get_num_vocab();

	for(int len : each_doc_len)
		total_words += len;

}

void Trainer::read_doc(string &file_path){
	// Add Document
	int doc_id = W.size();
	//cout << "     Doc id: " << doc_id << endl;
	W.push_back(vector<size_t>());
	each_doc_len.push_back(0); // doc_len initialize

	// Link filename and doc_id
	vector<string> components;
	split_string_by(file_path, '/', components);	// split by "/"
	string filename = components.back();
	doc_filename_to_id[filename] = doc_id;
	doc_id_to_filename[doc_id] = filename;

	// Read File
	ifstream ifs(file_path.c_str());
	string sentence;
	vector<string> sentences;
	while (getline(ifs, sentence) && !sentence.empty()){
		sentences.push_back(sentence);
	}
	for(string &sentence: sentences){
		vector<string> words;
		split_word_by(sentence, L' ', words);	// split by space
		add_sentence_to_doc(words, doc_id);
	}
}

void Trainer::add_sentence_to_doc(vector<string> &words, int doc_id){
	each_doc_len[doc_id] += words.size();

	for(string word: words){
		W[doc_id].push_back( vocab_c -> get_wordid(word) );
	}
}

////////////////
// Read Seeds //
////////////////
void Trainer::read_seed(string &seed_path){
	// Open file
	ifstream ifs(seed_path.c_str());
	string each_line;
	vector<string> lines;
	while(getline(ifs, each_line) && !each_line.empty()){
		lines.push_back(each_line);
	}

	// Read Lines
	num_topics = 0;
	for(string &line : lines){
		++num_topics; // # of seed topic should be equal to # of regular topic

		vector<string> words;
		split_word_by(line, L' ', words);

		// Read each seed
		int count = 0;
		unordered_map<size_t, double> phi_sk;
		for(string word : words){
			size_t word_id = vocab_c -> get_wordid(word);
			phi_sk[word_id] = 0.0;
			seed_words.insert(word_id);
			++count;
			++num_seed_vocab;
		}

		//// Set probabilities
		double prob = 1.0 / count;
		for(auto &p : phi_sk){
			p.second = prob;
		}

		phi_s.push_back(phi_sk);
		seed_num.push_back(count);
	}
}

// Summary Output
void Trainer::read_data_summary(){
	Rcout << "Data Summary:\n" << "  Number of documents: " << num_doc <<
		"\n  Number of all words in documents: " << total_words <<
		"\n  Number of unique words: " << num_vocab <<
		"\n  Number of unique seed words: " << num_seed_vocab <<
		"\n  Number of topics: " << num_topics << endl;
}

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
	slice_sampling_alpha(1e-9, 20.0, 3000);
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


List Trainer::get_output(int &show_words_num, bool &full_output){

	NumericVector iter_rvec;
	for(int iter : iter_vec)
		iter_rvec.push_back(iter);

	NumericVector loglik_rvec;
	for(double loglik : loglik_vec)
		loglik_rvec.push_back(loglik);

	NumericVector perplexity_rvec;
	for(double perplexity : perplexity_vec)
		perplexity_rvec.push_back(perplexity);

	// Z
	List Z_rlist;
	make_Z_rlist(Z_rlist);

	// X
	List X_rlist;
	make_X_rlist(X_rlist);

	// Collect Top ID
	List TopWords;
	make_TopWords(TopWords, show_words_num);

	// Iteration Info (alpha, n_x0, etc)
	List IterInfo;
	make_iter_info(IterInfo);

	if(!full_output){
		List return_list = List::create(Named("iter") = iter_rvec,
				Named("loglik") = loglik_rvec,
				Named("perplexity") = perplexity_rvec,
				Named("Z") = Z_rlist,
				Named("X") = X_rlist,
				Named("TopWords") = TopWords,
				Named("IterInfo") = IterInfo);

		return return_list;
	}

	// Create Additional Information
	// Word and Word ID
	List RawWords;
	List WordIDs;
	make_Words(RawWords, WordIDs);


	// Create Return Object
	List return_list = List::create(Named("iter") = iter_rvec,
			Named("loglik") = loglik_rvec,
			Named("perplexity") = perplexity_rvec,
			Named("Z") = Z_rlist,
			Named("X") = X_rlist,
			Named("TopWords") = TopWords,
			Named("RawWords") = RawWords,
			Named("WordIDs") = WordIDs,
			Named("IterInfo") = IterInfo);

	return return_list;
}

void Trainer::make_iter_info(List &IterInfo){
	IterInfo = List::create(Named("alpha") = iter_log::iter_alpha,
												  Named("p_k") = iter_log::iter_p_k,
													Named("n_x0_k") = iter_log::iter_n_x0_k,
													Named("n_x1_k") = iter_log::iter_n_x0_k);
}

void Trainer::make_Z_rlist(List &Z_rlist){
	for(vector<int> doc_z : Z){
		NumericVector docZ_r;
		for(int z : doc_z){
			docZ_r.push_back(z);
		}
		Z_rlist.push_back(docZ_r);
	}
}
void Trainer::make_X_rlist(List &X_rlist){
	for(vector<int> doc_x : X){
		NumericVector docX_r;
		for(int x : doc_x){
			docX_r.push_back(x);
		}
		X_rlist.push_back(docX_r);
	}
}
void Trainer::make_TopWords(List &TopWords, int &show_words_num){
	vector< vector< pair<string,int> > > wordscount; // words by topic
	make_wordscount(wordscount);
	sort_wordscount(wordscount);
	make_TopWordsRList(wordscount, TopWords, show_words_num);

}
void Trainer::make_wordscount(vector< vector< pair<string,int> > > &wordscount){
	size_t word_id;
	int z;
	string word;
	vector< pair<string,int> >::iterator ite;

	for(int k=0; k<num_topics; ++k)
		wordscount.push_back(vector< pair<string, int> >());

	for(int doc_id=0; doc_id<each_doc_len.size(); ++doc_id){
		//int doc_id_int = static_cast<int>(doc_id);

		for(int word_pos=0; word_pos<each_doc_len[doc_id]; ++word_pos){
			word_id =  W[doc_id][word_pos];
			word = vocab_c -> get_word(word_id);

			z = Z[doc_id][word_pos];

			ite = find_if(wordscount[z].begin(), wordscount[z].end(), comp(word));
			if (ite != wordscount[z].end()){
				int &c = ite->second;
				++c; // add count
			}else{
				// appear for the first time
				wordscount[z].push_back( make_pair(word, 1) );
			}
		}// for word_pos
	} // for doc_id
}
void Trainer::sort_wordscount(vector< vector< pair<string,int> > > &wordscount){
	for(int k=0; k<num_topics; ++k){
		// Sort by the second element, descending
		sort(wordscount[k].begin(), wordscount[k].end(), sortbysec_descend);
	}
}
void Trainer::make_TopWordsRList(vector< vector< pair<string,int> > > &wordscount,
		List &TopWords, int &show_words_num)
{
	double prop;

	for(int k=0; k<num_topics; ++k){
		IntegerVector words_count_r;
		NumericVector words_prop_r;
		CharacterVector words_r;

		for(int show=0; show<show_words_num; ++show){
			words_r.push_back( wordscount[k][show].first );
			words_count_r.push_back( wordscount[k][show].second );

			prop = (double)wordscount[k][show].second / (double)n_dk.col(k).sum();
			words_prop_r.push_back(prop);
		}
		List temp = List::create(Named("words") = words_r,
														 Named("count") = words_count_r,
														 Named("prop") = words_prop_r);
		TopWords.push_back(temp);
	}
}

void Trainer::make_Words(List &RawWords, List &WordIDs){
	int int_wid;
	string word;

	for(vector<size_t> doc_w : W){
		CharacterVector docW_r;
		NumericVector docID_r;

		for(size_t wid : doc_w){
			int_wid = static_cast<int>(wid);
			word = vocab_c -> get_word(wid);
			docID_r.push_back(int_wid);
			docW_r.push_back(word);
		}

		RawWords.push_back(docW_r);
		WordIDs.push_back(docID_r);
	}
}

// slice sampling

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




///////////////
// Read Data //
///////////////
// read_data.h


////////////////////
// Initialization //
////////////////////
// initialization.h

///////////////
// Iteration //
///////////////
// iteration.h

////////////////////////////////
// Calc_loglik and perplexity //
////////////////////////////////
// loglik.cpp

/////////////////
// Make output //
/////////////////
// output.cpp
