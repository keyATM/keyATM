#pragma once
using namespace std;
#include "trainer.hpp"

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
	vector< pair<string,int>>::iterator ite;

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

