#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
using namespace std;

struct faread {
	string hdr;
	string seq;
	int pair;
};
/*
void faread::operator=(faread r1, faread r2){
	r2.hdr = r1.hdr
	r2.seq = r1.seq;
	r2.qual = r1.qual;
}
*/
int main (int argc, char** argv) {
	
	if (argc == 1) {
		cerr << "Usage: " << argv[0] << "  <in.fasta>" << endl;
		return 0;
	}
	
	string inFile = argv[1];	
	string base = inFile.substr(0, int(inFile.find(".fa")));	
	faread read;
	
	ifstream in(inFile.c_str());
	string hdr;
	string seq;
	map<string,faread> pair1;
	map<string,faread> pair2;
	string key;
	while (in.good()){
		in >> hdr;
		in >> seq;
		hdr = hdr.substr(1);
		
		read.hdr = hdr;
		read.seq = seq;
		key = hdr.substr(0,hdr.length()-1);
		if (hdr.at(hdr.length()-1) == '1') {
			pair1[key] = read;
		} else {
			pair2[key] = read;
		}
	}

	ofstream p1out((base+"_p1.fasta").c_str());
	ofstream p2out((base+"_p2.fasta").c_str());
	ofstream upout((base+"_up.fasta").c_str());
	
	map<string,faread>::iterator it;
	faread*  tmp;
	for (it=pair1.begin(); it!=pair1.end(); it++) {
		if (pair2.find(it->first) != pair2.end()) {
			tmp =  & it->second;
			p1out << ">" << tmp->hdr << "\n" << tmp->seq << endl;
			tmp = & pair2.find(it->first)->second;
			p2out << ">" << tmp->hdr << "\n" << tmp->seq << endl;
			pair2.erase(it->first);
		} else {
			tmp = & it->second;
			upout << ">" << tmp->hdr << "\n" << tmp->seq << endl;
		}
	}
	
	for (it=pair2.begin(); it!=pair2.end(); it++) {
		tmp = & it->second;
		upout << ">" << tmp->hdr << "\n" << tmp->seq << endl;
	}
	
	p1out.close();
	p2out.close();
	upout.close();
	
}
