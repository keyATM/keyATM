#pragma once
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
using namespace Eigen;


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

int Bern(double &prob0, double &prob1){
	return randgen::bernoulli(prob1);
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

// Iteration info

// Others
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

vector<int> randomized_id_vec(int &num){
	return randgen::randomized_id_vec(num);
}
