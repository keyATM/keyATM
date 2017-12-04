#pragma once
#include <unordered_map>
#include <string>
#include <fstream>
using namespace std;

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
