#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <deque>
#include <map>
using namespace std;

int WIN_LEN = 21;
double wt;

map<char,int> tot;
int A = 0;
int T = 0;
int G = 0;
int C = 0;

double curr_gc;

deque<char> window; 

void printWindow(deque<char>& win, ostream& out){
	for(deque<char>::iterator it = win.begin(); it != win.end(); it++){
		out << *it;
	} 
	out << '\n';
}

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
		case 'g': G++; tot['g']++; break;
		case 'G': G++; tot['g']++; break;
		case 'c': C++; tot['c']++; break;
		case 'C': C++; tot['c']++; break;
		default: cerr << "Unrecognizable character: " << b << endl;
	}
}

void rembase(char& b) {
	switch (b) {
		case 'a': A--; break;
		case 'A': A--; break;
		case 't': T--; break;
		case 'T': T--; break;
		case 'g': G--; break;
		case 'G': G--; break;
		case 'c': C--; break;
		case 'C': C--; break;
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


void updateGC(char& o, char& n) {
	addbase(n);
	rembase(o);
	curr_gc = ((double)C+G)/((double)(WIN_LEN));
}

void finishSeq(deque<char>& win, list<double>& prof){
	double perc;
	for (int i = 0; i < WIN_LEN/2; i++){
		rembase(win.front());
		win.pop_front();
		perc = C+G;
		perc = perc/win.size();
		prof.push_back(perc);
	}
	while (win.size()>0 )
		win.pop_front();
	resetCounts();
}

void initlzSeq(deque<char>& win, list<double>& prof, ifstream& in){
	char b;
	for (int i = 0; i <= WIN_LEN/2; i++){
		b = (char) in.get();
		while (!isBase(b)){
			cerr << "Found non-base character where I'm not supposed to: >" << b << "<"<< endl;
			b = in.get();
		}	
		addbase(b);
		win.push_back(b);	
	}
	double perc = C+G;
	perc = perc/win.size(); 
	prof.push_back(perc); // push GC content for 1st position
	for (int i = 0; i < WIN_LEN/2; i++){
		b = (char) in.get();
		addbase(b);
		win.push_back(b);
		perc = C+G;
		perc = perc/win.size();
		prof.push_back(perc);
	}	
	curr_gc = perc;
} 
int main (int argc, char** argv) {
	
	if (argc == 1) {
		cerr << "Usage: gc_calc <in.fasta> <window_size> > <out.gc>" << '\n'
			 << "window_size should be odd. If even, window_size + 1 will be used." << endl;
		return 0;
	}	
		
	tot['a']=0;
	tot['t']=0;
	tot['g']=0;
	tot['c']=0;

	string inFile = argv[1];
	if (argc == 3) {
		WIN_LEN = atoi(argv[2]);
		if (WIN_LEN % 2 == 0) {
			WIN_LEN++;
		}
	}
	ifstream in(inFile.c_str());
	wt = 1.0;
	wt = 1.0/WIN_LEN;
	
	list<string> seq_names;
	list< list<double> > profiles;
	list<double> * tmp_prof = new list<double>();
	string hdr;
	char str[2048]; 
	const char* split = " \t|:";
	in.get();
	in.getline(str,2048,'\n');
	hdr = str;
	size_t spc_idx = hdr.find_first_of(split);
	hdr = hdr.substr(0,spc_idx);
	cerr << "counting sequence " << hdr << " spc_idx = " << spc_idx << endl;
	initlzSeq(window,*tmp_prof,in);	
	for(char b = in.get(); in.good(); b = in.get()){
		if (b == '>' ){ 
			finishSeq(window,*tmp_prof);
			seq_names.push_back(hdr);
			profiles.push_back(*tmp_prof);
			in.getline(str,2048,'\n');
			hdr = str;
			spc_idx = hdr.find_first_of(split);
			hdr = hdr.substr(0,spc_idx);
			cerr << "counting sequence " << hdr << " spc_idx = " << spc_idx << endl;
			tmp_prof = new list<double>();
			initlzSeq(window,*tmp_prof,in);
			continue;
		} else if (b == '\n')
			continue;
		updateGC(window.front(),b);
		window.pop_front();
		window.push_back(b);
		tmp_prof->push_back(curr_gc);
	}
	finishSeq(window,*tmp_prof);
	seq_names.push_back(hdr);
	profiles.push_back(*tmp_prof);
	
	list< list<double> >::iterator prof_it = profiles.begin();
	list<string>::iterator name_it = seq_names.begin(); 
	list<double>::iterator base_it;
	int pos = 1;
	while(prof_it != profiles.end()){
		for(base_it = prof_it->begin(); base_it != prof_it->end(); base_it++) {
			cout << *name_it <<'\t'<< pos <<'\t'<< *base_it << endl;
			pos++;
		}
		prof_it++;
		name_it++;
		pos = 1;
	}	
	double tot_gc = tot['g']+tot['c'];
	int tot_bases = tot['g']+tot['c']+tot['a']+tot['t'];
	cerr << "Total bases:" <<tot_bases<< endl;
	cerr << "   A: " << tot['a'] << endl; 
	cerr << "   T: " << tot['t'] << endl; 
	cerr << "   G: " << tot['g'] << endl; 
	cerr << "   C: " << tot['c'] << endl; 
	tot_gc = 100*tot_gc/tot_bases;	
	fprintf(stderr,"Genome-wide GC-content: %.2g%%\n",tot_gc);
}

