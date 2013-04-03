#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <deque>
#include <map>
#include <cmath>
using namespace std;

double wt;
int num_amb;

map<char,int> tot;
int A = 0;
int T = 0;
int G = 0;
int C = 0;

double win_G = 0;
double win_C = 0;
double win_size = 0;

double get_and_reset_window() {
	double ret = win_G+win_C;
	ret = ret/win_size;
	win_G = 0;
	win_C = 0;
	win_size = 0;
	return ret;
}

double curr_gc;
void resetCounts(){
	A = 0;
	T = 0;
	G = 0;
	C = 0;
	curr_gc = 0.0;
}

bool isGC(char& b){
	return b=='c' || b=='C' || b=='g'|| b == 'G';
}


void addbase(char& b) {
	switch (b) {
		case 'a': A++; tot['a']++; break;
		case 'A': A++; tot['a']++; break;
		case 't': T++; tot['t']++; break;
		case 'T': T++; tot['t']++; break;
		case 'g': win_G++; G++; tot['g']++; break;
		case 'G': win_G++; G++; tot['g']++; break;
		case 'c': win_C++; C++; tot['c']++; break;
		case 'C': win_C++; C++; tot['c']++; break;
		case 'M': num_amb++; break;
		case 'm': num_amb++; break;
		case 'R': num_amb++; break;
		case 'r': num_amb++; break;
		case 'W': num_amb++; break;
		case 'w': num_amb++; break;
		case 'S': num_amb++; break;
		case 's': num_amb++; break;
		case 'Y': num_amb++; break;
		case 'y': num_amb++; break;
		case 'K': num_amb++; break;
		case 'k': num_amb++; break;
		case 'V': num_amb++; break;
		case 'v': num_amb++; break;
		case 'H': num_amb++; break;
		case 'h': num_amb++; break;
		case 'D': num_amb++; break;
		case 'd': num_amb++; break;
		case 'B': num_amb++; break;
		case 'b': num_amb++; break;
		case 'X': num_amb++; break;
		case 'x': num_amb++; break;
		case 'n': num_amb++; break;
		case 'N': num_amb++; break;
		default: cerr << "Unrecognizable character: " << b << endl;
	}
}

bool isBase(char& b) {
	
	switch (b) {
		case 'a': return true;
		case 'A': return true;
		case 't': return true;
		case 'T': return true;
		case 'g': return true;
		case 'G': return true;
		case 'c': return true;
		case 'C': return true;
		default: return false; 
	}
}

int main (int argc, char** argv) {
	
	if (argc == 1) {
		cerr << "Usage: calc_avg-gc <in.fasta> > <out.gc>" << '\n';
		return 0;
	}	
		
	tot['a']=0;
	tot['t']=0;
	tot['g']=0;
	tot['c']=0;
	num_amb = 0;
	string inFile = argv[1];
	ifstream in(inFile.c_str());
	list<string> seq_names;
	list<double> gc_cont;
	list<double> window_gc;
	string hdr;
	const char* split = " \t|:";
	in.get();
	getline(in,hdr);
	size_t spc_idx = hdr.find_first_of(split);
	hdr = hdr.substr(0,spc_idx);
//	cerr << "counting sequence " << hdr << " spc_idx = " << spc_idx << endl;
	int len = 0;
	double GC = 0;
	for(char b = in.get(); in.good(); b = in.get()){
		if (b == '>' ){ 
//			cerr << hdr << " -> G: "<< G << "  C:  " << C << "  Length: " << len << endl;
			GC = G+C;
			GC = GC/len;
			seq_names.push_back(hdr);
			gc_cont.push_back(GC);
			resetCounts();
			len = 0;
			getline(in,hdr);
			spc_idx = hdr.find_first_of(split);
			hdr = hdr.substr(0,spc_idx);
			continue;
		} else if (b == '\n')
			continue;
		addbase(b);
		len++;
		win_size++;
		if (win_size == 1000) {
			window_gc.push_back(get_and_reset_window());
		}
	}
//	cerr << hdr << " -> G: "<< G << "  C:  " << C << "  Length: " << len << endl;
	GC = G+C;
	GC = GC/len;
	seq_names.push_back(hdr);
	gc_cont.push_back(GC);
	
	list<string>::iterator ctg_it = seq_names.begin();
	list<double>::iterator gc_it = gc_cont.begin();
	
	while(ctg_it != seq_names.end()){
		cout << *ctg_it << '\t' << *gc_it << endl;
		ctg_it++;
		gc_it++;
	}
	
	list<double>::iterator win_it = window_gc.begin();
	double s1 = 0.0;
	double s2 = 0.0;
	while (win_it != window_gc.end()){
		s1 += *win_it;
		s2 += pow(*win_it,2);
		win_it++;
	}
	double stdev = sqrt((s2 - ((pow(s1,2))/window_gc.size()))/(window_gc.size()-1)) * 100;
	double mean = s1/window_gc.size() * 101;

	double tot_gc = tot['g']+tot['c'];
	int tot_bases = tot['g']+tot['c']+tot['a']+tot['t'];
	cerr << "Total bases:" <<tot_bases<< endl;
	cerr << "   A: " << tot['a'] << endl; 
	cerr << "   T: " << tot['t'] << endl; 
	cerr << "   G: " << tot['g'] << endl; 
	cerr << "   C: " << tot['c'] << endl; 
	tot_gc = 100*tot_gc/tot_bases;	
	fprintf(stderr,"Genome-wide GC-content: %.2f +/- %.2f%%\n",tot_gc,stdev);
	if (num_amb > 0) {
		cerr << "Ignored " << num_amb << " ambiguous bases." << endl;
	}
}

