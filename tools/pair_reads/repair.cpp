#include <iostream>
#include <fstream>
#include <istream>
#include <map>
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace std;

bool fastq;
int pair_reads_call = 0;
double p = 1.0;
map<string,struct read> pair1;
map<string,struct read> pair2;

struct read {
	string hdr;
	string seq;
	string qual;
	int pair;
};


struct printer {
	ofstream* p1; 
	ofstream* p2; 
	bool print_fasta;
	bool revcomp;
		
	printer(string prefix, string base, string suffix, bool print_fasta, bool revcomp, bool shuf){
		if (shuf) {
			p1 = new ofstream((prefix+base+"_shuf"+suffix).c_str());
			p2 = p1; 
		} else {
			p1 = new ofstream((prefix+base+"_p1"+suffix).c_str());
			p2 = new ofstream((prefix+base+"_p2"+suffix).c_str());
		}
	}
	
	void print(read& r1, read& r2){
		if (fastq && !print_fasta) {
			p1out << "@" << r1->hdr << "\n" << r1->seq << "\n+" << r1->hdr << "\n" << r1->qual << endl;
			p2out << "@" << r2->hdr << "\n" << r2->seq << "\n+" << r2->hdr << "\n" << r2->qual << endl;
		} else {
			p1out << ">" << r1->hdr << "\n" << r1->seq << endl;
			p2out << ">" << r2->hdr << "\n" << r2->seq << endl;

		}
	
	}

}


char comp(char b) {
    switch(b){
        case 'A': return 'T';
        case 'a': return 't';
        case 'T': return 'A';
        case 't': return 'a';
        case 'G': return 'C';
        case 'g': return 'c';
        case 'C': return 'G';
        case 'c': return 'g';
        case 'M': return 'K';
        case 'm': return 'k';
        case 'R': return 'Y';
        case 'r': return 'y';
        case 'W': return 'W';
        case 'w': return 'w';
        case 'S': return 'S';
        case 's': return 's';
        case 'Y': return 'R';
        case 'y': return 'r';
        case 'K': return 'M';
        case 'k': return 'm';
        case 'V': return 'B';
        case 'v': return 'b';
        case 'H': return 'D';
        case 'h': return 'd';
        case 'D': return 'H';
        case 'd': return 'h';
        case 'B': return 'V';
        case 'b': return 'v';
        case 'X': return 'X';
        case 'x': return 'x';
        case 'N': return 'N';
        case 'n': return 'n';
        default : { cerr << "Unrecognized character: " << b << endl;
                    return '-';}
    }
}

void rev_comp(struct read& r){
    int j = -1;
    int stop = -1;
    if (r.seq.length() % 2 == 1)
        stop = r.seq.length()/2;
    else
        stop = r.seq.length()/2 - 1;

    for (int i = 0; i < stop; i++){
        j = r.seq.length()-1-i;
        char to_j = comp(r.seq.at(i));
        char to_i = comp(r.seq.at(j));
        r.seq.replace(i,1,1,to_i);
        r.seq.replace(j,1,1,to_j);
    }

    if (fastq) {
        for (int i = 0; i < stop; i++){
            j = r.qual.length()-1-i;
            char to_j = r.qual.at(i);
            char to_i = r.qual.at(j);
            r.qual.replace(i,1,1,to_i);
            r.qual.replace(j,1,1,to_j);
        }

    }
}


double unif_rand() {
	return rand()/double(RAND_MAX);
}

void pair_reads(istream& in, bool fastq){ 
	pair_reads_call++;
	struct read r;
	string hdr;
	string q_hdr;
	string seq;
	string qual;
	string key;
	while (in >> hdr){
		in >> seq;
		if (fastq) {
			in >> q_hdr;
			in >> qual;
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
		} else if (hdr.at(hdr.length()-1) == '2' || hdr.at(hdr.length()-1) == '3') {
			r.pair = 2;
			pair2[key] = r;
		} else {
			if (pair1.find(key) == pair1.end()){
				r.pair = 1;
				pair1[key] = r;
			} else {
				r.pair = 2;
				pair2[key] = r;
			}
		}
	}
}

void print_paired(string prefix, string base, string suffix, bool print_fasta, bool revcomp){
	ofstream p1out((prefix+base+"_p1"+suffix).c_str());
	ofstream p2out((prefix+base+"_p2"+suffix).c_str());
	ofstream upout;
	
	map<string, struct read>::iterator it;
	struct read*  tmp1;
	struct read*  tmp2;
	double U = 0.0;
	for (it=pair1.begin(); it!=pair1.end(); it++) {
		if (p != 1.0) {
			U = unif_rand();
			if (U > p) continue;
		}
		if (pair2.find(it->first) != pair2.end()) {
			tmp1 = & it->second;
			tmp2 = & pair2.find(it->first)->second; 
			if (revcomp) {
				rev_comp(*tmp1);
				rev_comp(*tmp2);
			}
			if (fastq && !print_fasta) {
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
			if (revcomp) {
				rev_comp(*tmp1);
			}
			if (fastq && !print_fasta)
				upout << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
			else{
				upout << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
			}
		}
	}
	
	for (it=pair2.begin(); it!=pair2.end(); it++) {
		if (p != 1.0) {
			U = unif_rand();
			if (U > p) continue;
		}
		if (!upout.is_open()){
			upout.open((prefix+base+"_up"+suffix).c_str());
		}
		tmp2 = & it->second;
		if (revcomp) {
			rev_comp(*tmp2);
		}
		if (fastq && !print_fasta)
			upout << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
		else {
			upout << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;
		}
	}
	
	p1out.close();
	p2out.close();
	if (upout.is_open())
		upout.close();
	
}

void print_shuffled(string prefix, string base, string suffix, bool print_fasta, bool revcomp){

	ofstream paired((prefix+base+"_shuf"+suffix).c_str());
	ofstream unpaired;
	
	map<string, struct read>::iterator it;
	struct read*  tmp1;
	struct read*  tmp2;
	double U = 0.0;	
	for (it=pair1.begin(); it!=pair1.end(); it++) {
		if (p != 1.0) {
			U = unif_rand();
			if (U > p) continue;
		}
		if (pair2.find(it->first) != pair2.end()) {
			tmp1 = & it->second;
			tmp2 = & pair2.find(it->first)->second; 
			if (revcomp){
				rev_comp(*tmp1);
				rev_comp(*tmp2);
			}
			if (fastq && !print_fasta) {
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
			if (revcomp){
				rev_comp(*tmp1);
			}
			if (fastq && !print_fasta)
				unpaired << "@" << tmp1->hdr << "\n" << tmp1->seq << "\n+" << tmp1->hdr << "\n" << tmp1->qual << endl;
			else {
				unpaired << ">" << tmp1->hdr << "\n" << tmp1->seq << endl;
			}
		}
	}
	
	for (it=pair2.begin(); it!=pair2.end(); it++) {
		if (p != 1.0) {
			U = unif_rand();
			if (U > p) continue;
		}
		if (!unpaired.is_open()) {
			unpaired.open((prefix+base+"_up"+suffix).c_str());
		}
		tmp2 = & it->second;
		if (revcomp){
			rev_comp(*tmp2);
		}
		if (fastq && !print_fasta)
			unpaired << "@" << tmp2->hdr << "\n" << tmp2->seq << "\n+" << tmp2->hdr << "\n" << tmp2->qual << endl;
		else { 
			unpaired << ">" << tmp2->hdr << "\n" << tmp2->seq << endl;
		}
	}
	
	paired.close();
	if (unpaired.is_open())
		unpaired.close();
	
}

void pipe_seq(istream& in, ostream& out, bool print_fasta, bool revcomp) {
	struct read r;
	if (in >> r.hdr){
		in >> r.seq;
		if (fastq) {
			in >> r.hdr;
			in >> r.qual;
			r.hdr.replace(0,1,1,'@');
		}
		if (revcomp) 
			rev_comp(r);

		if (fastq && !print_fasta){
			out << r.hdr << endl;
			out << r.seq << endl;
			r.hdr.replace(0,1,1,'+');
			out << r.hdr << endl;
			out << r.qual << endl;
		} else {
			r.hdr.replace(0,1,1,'>');
			out << r.hdr << endl;
			out << r.seq << endl;
		}
	}	
}

void skip_seq(istream& in, bool fastq) {
	string hdr;
	string seq;
	if (in >> hdr){
		in >> seq;
		if (fastq) {
			in >> hdr;
			in >> seq;
		}
	}
}

void sample_paired(ifstream& in1, ifstream& in2, ofstream& p1out, ofstream& p2out, bool fastq, bool print_fasta) {
	double U = 0.0;	
	while (in1.good() && in2.good()){
		if (p != 1.0) {
			U = unif_rand();
			if (U > p){
				skip_seq(in1,fastq);
				skip_seq(in2,fastq);
				continue;
			}	
		}
		pipe_seq(in1,p1out,fastq,print_fasta);
		pipe_seq(in2,p2out,fastq,print_fasta);
	}
}

void sample_shuffled(istream& in, ofstream& out, bool fastq, bool print_fasta){
	double U = 0.0;	
	while (in.good()){
		if (p != 1.0) {
			U = unif_rand();
			if (U > p){
				skip_seq(in,fastq);
				skip_seq(in,fastq);
				continue;
			}	
		}
		pipe_seq(in,out,fastq,print_fasta);
		pipe_seq(in,out,fastq,print_fasta);
	}
}

void shuffle_paired(ifstream& in1, ifstream& in2, ofstream& out, bool fastq, bool print_fasta) {
	double U = 0.0;	
	while (in1.good() && in2.good()){
		if (p != 1.0) {
			U = unif_rand();
			if (U > p){
				skip_seq(in1,fastq);
				skip_seq(in2,fastq);
				continue;
			}	
		}
		pipe_seq(in1,out,fastq,print_fasta);
		pipe_seq(in2,out,fastq,print_fasta);
	}
}

void split_shuffled(istream& in, ofstream& p1out, ofstream& p2out, bool print_fasta, bool revcomp){
	double U = 0.0;	
	while (in.good()){
		if (p != 1.0) {
			U = unif_rand();
			if (U > p){
				skip_seq(in,fastq);
				skip_seq(in,fastq);
				continue;
			}	
		}
		pipe_seq(in,p1out,fastq,print_fasta);
		pipe_seq(in,p2out,fastq,print_fasta);
	}
}


void usage(const char* name){
	cout << "Usage: "<<name<<" [options] <base> <reads.in>"<<endl;
	cout << " where "<< endl;
	cout << "       <base>       basename for output files\n";
	cout << "       <reads.in>   a list of files containing all reads to re-pair\n";
	cout << "                    Files can be in fastq or fasta format, but all\n"; 
	cout << "                    files must be in the same format. If no read files \n";
	cout << "                    are listed, input is read from standard in.\n"; 
	cout << " options:\n";
	cout << "        -p <string>    a prefix to add to <base>. This can be used to\n"; 
	cout << "                       specify an output directory.\n";
	cout << "        -s <string>    the suffix to append to the output files.\n"; 
	cout << "        --paired       assume reads are paired.\n";
	cout << "        --shuf         print pairs in one file where paired reads are printed\n";
	cout << "                       on consecutive lines (a.k.a. shuffled).\n";
	cout << "        --split        split pairs from one or more shuffled files into two files\n";
	cout << "                       Note: this option is only applicable when reads are already paired.\n";
	cout << "        -r <double>    randomly sample reads or read-pairs with specified probability\n"; 
	cout << "        --seed <long>  seed for randomly sampling reads\n";	
	cout << "        --rev          reverse complement each sequence.\n";
	cout << "        --fasta        output in fasta format.\n";
	cout << "        --quiet        do not print progress messages.\n";
	cout << "        --debug        run in debug mode.\n\n";
	 
}


int main (int argc, const char** argv) {
	if (argc == 1) {
		usage(argv[0]);
		return 0;
	}
	string prefix = "";
	string suffix = "";
	string base = "";
	unsigned int seed = time(NULL);
	//bool fastq;
	bool revcomp = false;
	bool output_fasta=false;
	bool shuffle = false;
	bool split = false;
	bool paired = false;
	bool debug = false;
	bool quiet = false;
	int start = 1;
	int i = 1;
	while (argv[i][0] == '-') {
		if (argv[i][1]=='s'){
			suffix = argv[++i];
			start+=2;
		} else if (argv[i][1]=='p') {
			prefix = argv[++i];
			start+=2;
		} else if (argv[i][1]=='r') {
			p = atof(argv[++i]);
			start+=2;
		} else if (argv[i][1]=='-'){
			if (strcmp(argv[i],"--shuf")==0){
				shuffle = true;
			} else if (strcmp(argv[i],"--paired")==0){
				paired = true;
			} else if (strcmp(argv[i],"--split")==0){
				split = true;
			} else if (strcmp(argv[i],"--quiet")==0){
				quiet = true;
			} else if (strcmp(argv[i],"--debug")==0){
				debug = true;
			} else if (strcmp(argv[i],"--fasta")==0){
				output_fasta = true;
			} else if (strcmp(argv[i],"--seed")==0){
				seed = atoi(argv[++i]);
				start++;
			}
			start++;
		} else {
			cerr << "Unrecognized argument: " << argv[i] << endl;
		}
		i++;
	}
	srand(seed);
	if (p != 1.0) cout << "subsampling reads with seed " << seed << endl;
	if (debug) cerr << argc << " arguments. Starting at " << start << endl;
	
	char c;
	istream* in;
	base = argv[start++];	
	if (fopen(base.c_str(),"r")) {
		cerr << "Missing <base> argument: " << base << " is a file.\n";
		usage(argv[0]);
		return 0;
	}
	if (paired) { 
		if (shuffle) { // need two files
			if (argc - start == 2) {
				ifstream in1(argv[start++]);
				ifstream in2(argv[start++]);
				ofstream out((prefix+base+"_shuf"+suffix).c_str());
				fastq = in1.peek()=='@';
				if (!quiet) cout << "Shuffling paired reads from " << argv[start-2] << " and " << argv[start-1] << "\n";
				shuffle_paired(in1,in2,out,output_fasta,revcomp);
			} else {
				cerr << "Two files are required for shuffling paired reads.\n";
				usage(argv[0]);
			}
		} else if (split) {  
			if (!quiet) cout << "Splitting shuffled reads from";
			ofstream p1out((prefix+base+"_p1"+suffix).c_str());
			ofstream p2out((prefix+base+"_p2"+suffix).c_str());
			if (argc - start == 0) {
				in = &cin;
				fastq = in->peek() == '@';
				if (!quiet) cout << " standard input.\n";
				split_shuffled(*in,p1out,p2out,output_fasta,revcomp);
			} else {
				do {			
					if (!quiet) cout << "\n\t" << argv[start];
					in = new ifstream(argv[start],ifstream::in);
					fastq = in->peek() == '@';
				//	if (debug) cerr << start << "  " << argv[start] << endl;
					split_shuffled(*in,p1out,p2out,output_fasta,revcomp);
					delete in;
					start++;
				} while (start < argc);
				if (!quiet) cout << '\n';
				//	fb.open(argv[i],ios::in);
				//	in = new istream(&fb);
				//	split_shuffled(*in,p1out,p2out);
				
			}
			return 1;
		} else if (p < 1.0) { // just randomly sample reads 
			if (argc - start == 2) { // assume two files are in non-shuffled paired format
				ofstream p1out((prefix+base+"_p1"+suffix).c_str());
				ofstream p2out((prefix+base+"_p2"+suffix).c_str());
				ifstream in1(argv[start++]);
				ifstream in2(argv[start++]);
				sample_paired(in1,in2,p1out,p2out,fastq,output_fasta);
			} else { // assume all files are shuffled
				ofstream out ((prefix+base+"_shuf"+suffix).c_str());
				ifstream* in; 
				while (start < argc) {
					in = new ifstream(argv[start++]);
					sample_shuffled(*in,out,fastq,output_fasta);
					delete in;
				}
			}

		} else {
			if (!quiet) cout << "Doing nothing\n";
		}
	} else  { // reads need to be re-paired using a hash
		if (!quiet) cout << "Pairing reads from" ;
		if (argc - start == 0) {
			in = &cin;
			fastq = in->peek() == '@';
			if (!quiet) cout << " standard input.\n";
			pair_reads(*in,fastq);
		} else {
			do {
				if (!quiet) cout << "\n\t" << argv[start-1];
				in = new ifstream(argv[start++],ifstream::in);
				fastq = in->peek() == '@';
				pair_reads(*in,fastq);
				delete in;
			} while (start < argc);
			if (!quiet) cout << '\n';
//			for (int i = start; i < argc; i++){
//				cerr << "Loading reads from " << argv[i] << endl;
//				fb.open(argv[i],ios::in);
//				in = new istream(&fb);
//				pair_reads(*in);
//				fb.close();
//				delete in;
//			}
		}
		if (shuffle)
			print_shuffled(prefix,base,suffix,output_fasta,revcomp);
		else
			print_paired(prefix,base,suffix,output_fasta,revcomp);
	}

}
