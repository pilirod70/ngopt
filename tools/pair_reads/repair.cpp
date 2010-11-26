#include <iostream>
#include <fstream>
#include <istream>
#include <map>
#include <cstdio>
#include <cstring>

using namespace std;

struct read {
	string hdr;
	string seq;
	string qual;
	int pair;
};


bool fastq;
int load_reads_call = 0;
map<string,struct read> pair1;
map<string,struct read> pair2;

void load_reads(istream& in){ 
	load_reads_call++;
	struct read r;
	string hdr;
	string q_hdr;
	string seq;
	string qual;
	string key;
	while (in.good()){
		in >> hdr;
		in >> seq;
		if (fastq) {
			in >> q_hdr;
			in >> qual;
		}
		if (!in.good()) {
			break;
		}
		hdr = hdr.substr(1);
		
		r.hdr = hdr;
		r.seq = seq;
		
		if (fastq){
			r.qual = qual;
		}
		key = hdr.substr(0,hdr.length()-2);

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
	ofstream upout;
	
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
			if (!upout.is_open()){
				upout.open((prefix+base+"_up"+suffix).c_str());
			}
			tmp1 = & it->second;
			if (fastq)
				upout << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
			else
				upout << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
		}
	}
	
	for (it=pair2.begin(); it!=pair2.end(); it++) {
		if (!upout.is_open()){
			upout.open((prefix+base+"_up"+suffix).c_str());
		}
		tmp2 = & it->second;
		if (fastq)
			upout << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
		else 
			upout << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;
	}
	
	p1out.close();
	p2out.close();
	if (upout.is_open())
		upout.close();
	
}

void print_shuffled(string prefix, string base, string suffix){

	ofstream paired((prefix+base+"_shuf"+suffix).c_str());
	ofstream unpaired;
	
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
			if (!unpaired.is_open()) {
				unpaired.open((prefix+base+"_up"+suffix).c_str());
			}
			tmp1 = & it->second;
			if (fastq)
				unpaired << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
			else
				unpaired << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
		}
	}
	
	for (it=pair2.begin(); it!=pair2.end(); it++) {
		if (!unpaired.is_open()) {
			unpaired.open((prefix+base+"_up"+suffix).c_str());
		}
		tmp2 = & it->second;
		if (fastq)
			unpaired << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
		else 
			unpaired << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;
	}
	
	paired.close();
	if (unpaired.is_open())
		unpaired.close();
	
}

void split_shuffled(istream& in, ofstream& p1out, ofstream& p2out){
	string buf;
	int stop = 2;
	if (fastq) 
		stop = 4;
	while (in.good()){
		for (int i = 0; i < stop; i++) {
			in >> buf;
			p1out << buf << endl;
		}
		for (int i = 0; i < stop; i++) {
			in >> buf;
			p2out << buf << endl;
		}
	}
	
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
	cout << "                    on consecutive lines (a.k.a. shuffled).\n";
	cout << "        --split     split a set of shuffled reads up into two files.\n\n";
	 
}


int main (int argc, char** argv) {
	cerr << argc << " arguments.\n";	
	if (argc == 1) {
		usage(argv[0]);
	//	cerr << "Usage: " << argv[0] << "  <reads.in>" << endl;
		return 0;
	}
	string prefix = "";
	string suffix = "";
	string base = "";
	bool shuffle = false;
	bool split = false;
	int start = 1;
	int i = 1;
	while (argv[i][0] == '-') {
		if (argv[i][1]=='s'){
			suffix = argv[++i];
			start+=2;
		} else if (argv[i][1]=='p') {
			prefix = argv[++i];
			start+=2;
		} else if (argv[i][1]=='-'){
		//	argv[i]+=2;
			if (strcmp(argv[i],"--shuf")==0){
				shuffle = true;
				cerr << "Shuffling reads" << endl;
			} else if (strcmp(argv[i],"--split")==0)
				split = true;
				cerr << "Splitting reads" << endl;
			start++;
		} else {
			cerr << "Unrecognized argument: " << argv[i] << endl;
		}
		i++;
	}
	char c;
	istream* in;
	base = argv[start++];	
	if (fopen(base.c_str(),"r")) {
		cerr << "Missing <base> argument\n";
		usage(argv[0]);
		return 0;
	}
	if (split) {
		cerr << "Splitting reads\n" ;
		ofstream p1out((prefix+base+"_p1"+suffix).c_str());
		ofstream p2out((prefix+base+"_p2"+suffix).c_str());
		if (argc - start == 0) {
			in = &cin;
			fastq = in->peek() == '@';
			cerr << "Loading reads from standard input." << endl;
			split_shuffled(*in,p1out,p2out);
		} else {
			filebuf fb;
			fb.open(argv[start++],ios::in);
			in = new istream(&fb);
			fastq = in->peek() == '@';
			cerr << start << "  " << argv[start] << endl;
			split_shuffled(*in,p1out,	p2out);
			cerr << "Loading reads from multiple files." << endl;
			for (int i = start; i < argc; i++){
				fb.open(argv[i],ios::in);
				in = new istream(&fb);
				split_shuffled(*in,p1out,p2out);
			}
		}
		return 1;
	} else{ // reads need to be re-paired using a hash
		cerr << "Pairing reads\n" ;
		if (argc - start == 0) {
			in = &cin;
			fastq = in->peek() == '@';
			cerr << "Loading reads from standard input." << endl;
			load_reads(*in);
		} else {
			filebuf fb;
			fb.open(argv[start++],ios::in);
			in = new istream(&fb);
			fastq = in->peek() == '@';
			cerr << "Loading reads from multiple files." << endl;
			load_reads(*in);
			fb.close();
			delete in;
			for (int i = start; i < argc; i++){
				cerr << "Loading reads from " << argv[i] << endl;
				fb.open(argv[i],ios::in);
				in = new istream(&fb);
				load_reads(*in);
				fb.close();
				delete in;
			}
		}
		if (shuffle)
			print_shuffled(prefix,base,suffix);
		else
			print_paired(prefix,base,suffix);
	}

}
