#include <cstdio>
#include <list>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <cstdlib>

using namespace std;

struct fa_entry {
	string name;
	list<char> seq;
};


void print_seq(struct fa_entry& fa, ostream& out) {
	out << '>' << fa.name << endl;
	list<char>::iterator it = fa.seq.begin(); 
	int len = 0;
	while (it != fa.seq.end()){
		out << *it;
		it++;
	}
	out << '\n';
}

bool is_unk(char c) {
	return c == '-' || c == 'N' || c == 'n' || c == 'X' || c == 'x';
}

bool is_nuc(char c) {
	return c == 'A' || c == 'T' || c == 'G' || c == 'C' || 
		   c == 'a' || c == 't' || c == 'g' || c == 'c' ; 
}

int main(int argc, char** argv) {
	
	if (argc != 3) {
		cout << "Usage: spscaf <mingap> <in.fasta>\n"
			 << "gaps longer than <mingap> will be split into two contigs\n"
			 << "Output printed to standard out\n";
		return -1;
	}

	int MINGAP = atoi(argv[1]);
	ifstream in(argv[2]);
	string infile(argv[2]);
	string base(infile.substr(0,infile.rfind(".fa")));	
	
	ofstream out_cis((base+".cis.txt").c_str());
	ofstream out_ctgs((base+".contigs.fasta").c_str());

	struct fa_entry* tmp = new fa_entry;
	char c = (int) in.get();
	if (c != '>') {
		cerr << "File not in fasta format\n";
		return -1;
	}
	cout << "Splitting scaffolds in " << argv[1] << endl
		 << "Gaps of length " << MINGAP << " or longer will be split\n";
	char* scaf;
 	char buf[256];	
	in.getline(buf,256,'\n');	
	scaf = buf;
	stringstream ss;
	int contig = 1;
	ss << "contig" << contig;
	tmp->name = ss.str();
	ss.str("");
	c = (char) in.get();
	while (in.good()){
		//if (!in.good()) break; 
		if (c == '>'){
			if (tmp->seq.size() > 1) {
				print_seq(*tmp,out_ctgs);
				out_cis << tmp->name << "\t" << buf << "\n";
			}
			contig++;
			delete tmp;
			tmp = new fa_entry;
			in.getline(scaf,256,'\n');	
		//	scaf = buf;
			ss << "contig" << contig;
			cout << contig << "\t";
			tmp->name = ss.str();
			cout << tmp->name << endl;
			ss.str("");
		} else if (is_unk(c)) {
			int num_n = 1;
			while ((is_unk(c) || c=='\n') && c != '>'){
				c = (char) in.get();
				if (is_unk(c)){
					num_n++;
				}	
			}
			if (num_n <= MINGAP) { 
				cerr << "Found small gap of length " << num_n << endl;
				for (int i = 0; i < num_n; i++){
					tmp->seq.push_back('N');
				}
				continue;
			}
			if (tmp->seq.size() > 1) {
				print_seq(*tmp,out_ctgs);
				out_cis << tmp->name << "\t" << buf << "\n";
			}
			contig++;
			delete tmp;
			tmp = new fa_entry;
			ss << "contig" << contig;
			cout << contig << "\t";
			tmp->name = ss.str();
			cout << tmp->name << endl;
			ss.str("");
			if (is_nuc(c)) {
				tmp->seq.push_back(c);
			}
		} else if (c == '\n') {
			c = (char) in.get();
			continue;
		}  else {
			tmp->seq.push_back(c);
		}
		c = (char) in.get();
	}
	
	if (tmp->seq.size() > 1) {
		print_seq(*tmp,out_ctgs);
		out_cis << tmp->name << "\t" << scaf << "\n";		
		cout << tmp->seq.size() << endl;
	}
	delete tmp;
}	
