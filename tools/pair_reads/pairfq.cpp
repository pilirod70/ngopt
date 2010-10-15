#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
using namespace std;

struct fqread {
	string hdr;
	string seq;
	string qual;
	int pair;
};
/*
void fqread::operator=(fqread r1, fqread r2){
	r2.hdr = r1.hdr
	r2.seq = r1.seq;
	r2.qual = r1.qual;
}
*/
int main (int argc, char** argv) {
	
	if (argc == 1) {
		cerr << "Usage: " << argv[0] << "  <in.fastq>" << endl;
		return 0;
	}
	
	string inFile = argv[1];	
	string base = inFile.substr(0, int(inFile.find(".fa")));	
	fqread read;
	
	ifstream in(inFile.c_str());
	string hdr;
	string seq;
	string qual;
	map<string,fqread> pair1;
	map<string,fqread> pair2;
	string key;
	while (in.good()){
		in >> hdr;
		in >> seq;
		in >> hdr;
		in >> qual;
		hdr = hdr.substr(1);
		
		read.hdr = hdr;
		read.seq = seq;
		read.qual = qual;
		key = hdr.substr(0,hdr.length()-1);
		if (hdr.at(hdr.length()-1) == '1') {
			pair1[key] = read;
		} else {
			pair2[key] = read;
		}
	}

	ofstream p1out((base+"_p1.fastq").c_str());
	ofstream p2out((base+"_p2.fastq").c_str());
	ofstream upout((base+"_up.fastq").c_str());
	
	map<string,fqread>::iterator it;
	fqread*  tmp;
	for (it=pair1.begin(); it!=pair1.end(); it++) {
		if (pair2.find(it->first) != pair2.end()) {
			tmp =  & it->second;
			p1out << "@" << tmp->hdr << "\n" << tmp->seq << "\n+" << tmp->hdr << "\n" << tmp->qual << endl;
			tmp = & pair2.find(it->first)->second;
			p2out << "@" << tmp->hdr << "\n" << tmp->seq << "\n+" << tmp->hdr << "\n" << tmp->qual << endl;
			pair2.erase(it->first);
		} else {
			tmp = & it->second;
			upout << "@" << tmp->hdr << "\n" << tmp->seq << "\n+" << tmp->hdr << "\n" << tmp->qual << endl;
		}
	}
	
	for (it=pair2.begin(); it!=pair2.end(); it++) {
		tmp = & it->second;
		upout << "@" << tmp->hdr << "\n" << tmp->seq << "\n+" << tmp->hdr << "\n" << tmp->qual << endl;
	}
	
	p1out.close();
	p2out.close();
	upout.close();
	
}
 










