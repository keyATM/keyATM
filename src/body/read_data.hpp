#pragma once
#include "trainer.hpp"
#include "corpus.hpp"

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


