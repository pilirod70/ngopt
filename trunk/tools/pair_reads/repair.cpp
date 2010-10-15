#include <iostream>
#include <fstream>
#include <istream>
#include <map>
#include <stdio.h>
using namespace std;

struct read {
	string hdr;
	string seq;
	string qual;
	int pair;
};


bool fastq;

map<string,struct read> pair1;
map<string,struct read> pair2;

void load_reads(istream& in){ 

	struct read r;
	string hdr;
	string seq;
	string qual;
	string key;
	while (in.good()){
		in >> hdr;
		in >> seq;
		if (fastq) {
			in >> hdr;
			in >> qual;
		}
		hdr = hdr.substr(1);
		
		r.hdr = hdr;
		r.seq = seq;
		if (fastq);
			r.qual = qual;
		key = hdr.substr(0,hdr.length()-1);
		if (hdr.at(hdr.length()-1) == '1') {
			r.pair = 1;
			pair1[key] = r;
		} else {
			r.pair = 2;
			pair2[key] = r;
		}
	}
}


void print_paired(string prefix, string base, string suffix){

	ofstream p1out((prefix+base+"_p1"+suffix).c_str());
	ofstream p2out((prefix+base+"_p2"+suffix).c_str());
	ofstream upout((prefix+base+"_up"+suffix).c_str());
	
	map<string, struct read>::iterator it;
	struct read*  tmp1;
	struct read*  tmp2;
	for (it=pair1.begin(); it!=pair1.end(); it++) {
		if (pair2.find(it->first) != pair2.end()) {
			tmp1 = & it->second;
			tmp2 = & pair2.find(it->first)->second; 
			if (fastq) {
				p1out << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
				p2out << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
			} else {
				p1out << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
				p2out << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;

			}
			pair2.erase(it->first);
		} else {
			tmp1 = & it->second;
			if (fastq)
				upout << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
			else
				upout << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
		}
	}
	
	for (it=pair2.begin(); it!=pair2.end(); it++) {
		tmp2 = & it->second;
		if (fastq)
			upout << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
		else 
			upout << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;
	}
	
	p1out.close();
	p2out.close();
	upout.close();
	
}

void print_shuffled(string prefix, string base, string suffix){

	ofstream paired((prefix+base+"_shuf"+suffix).c_str());
	ofstream unpaired((prefix+base+"_up"+suffix).c_str());
	
	map<string, struct read>::iterator it;
	struct read*  tmp1;
	struct read*  tmp2;
	for (it=pair1.begin(); it!=pair1.end(); it++) {
		if (pair2.find(it->first) != pair2.end()) {
			tmp1 = & it->second;
			tmp2 = & pair2.find(it->first)->second; 
			if (fastq) {
				paired << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
				paired << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
			} else {
				paired << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
				paired << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;

			}
			pair2.erase(it->first);
		} else {
			tmp1 = & it->second;
			if (fastq)
				unpaired << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
			else
				unpaired << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
		}
	}
	
	for (it=pair2.begin(); it!=pair2.end(); it++) {
		tmp2 = & it->second;
		if (fastq)
			unpaired << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
		else 
			unpaired << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;
	}
	
	paired.close();
	unpaired.close();
	
}

void usage(char* name){
	cout << "Usage: "<<name<<" [options] <base> <reads.in>"<<endl;
	cout << " where "<< endl;
	cout << "       <base>       basename for output files\n";
	cout << "       <reads.in>   a list of files containing all reads to re-pair\n";
	cout << "                    Files can be in fastq or fasta format, but all\n"; 
	cout << "                    files must be in the same format. If no read files \n";
	cout << "                    are listed, input is read from standard in.\n"; 
	cout << " options:\n";
	cout << "        -p <string> a prefix to add to <base>. This can be used to\n"; 
	cout << "                    specify an output directory.\n";
	cout << "        -s <string> the suffix to append to the output files.\n"; 
	cout << "        --shuf      print pairs in one file where paired reads are printed\n";
	cout << "                    on consecutive lines (a.k.a. shuffled).\n\n";
	 
}


int main (int argc, char** argv) {
	
	if (argc == 1) {
		usage(argv[0]);
	//	cerr << "Usage: " << argv[0] << "  <reads.in>" << endl;
		return 0;
	}
	string prefix = "";
	string suffix = "";
	string base = "";
	bool shuffle = false;
	int start = 1;
	int i = 1;
	while (argv[i][0] == '-') {
		if (argv[i][1]=='s'){
			i++;
			suffix = argv[i++];
			start+=2;
		} else if (argv[i][1]=='p') {
			i++;
			prefix = argv[i++];
			start+=2;
		} else if (argv[i][1]=='-'){
			argv[i]+=2;
			if (strcmp(argv[i],"shuf"))
				shuffle = true;
		} else {
			cerr << "Unrecognized argument: " << argv[i] << endl;
		}
	}
	char c;
	istream* in;
	base = argv[start++];
	if (argc - start == 0) {
		in = &cin;
		fastq = in->peek() == '@';
		load_reads(*in);
	} else {
		filebuf fb;
		fb.open(argv[start++],ios::in);
		in = new istream(&fb);
		fastq = in->peek() == '@';
		load_reads(*in);
		for (int i = start; i < argc; i++){
			fb.open(argv[i],ios::in);
			in = new istream(&fb);
			load_reads(*in);
		}
	}

	if (shuffle)
		print_shuffled(prefix,base,suffix);
	else
		print_paired(prefix,base,suffix);
}
