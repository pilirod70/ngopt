#include <cstdio>
#include <list>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>

using namespace std;

struct fa_entry {
	const char* name;
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
	return c == '-' || c =='\n' || c == 'N' || c == 'n' || c == 'X' || c == 'x';
}


int main(int argc, char** argv) {
	
	if (argc != 2) {
		cout << "Usage: spscaf <in.fasta>\n" << "Output printed to standard out\n";
		return -1;
	}

	ifstream in(argv[1]);
	string infile(argv[1]);
	string base(infile.substr(0,infile.find_last_of(".fa")));	
	
	ofstream out_cis((base+".contigs.fasta").c_str());
	ofstream out_ctgs((base+".cis.txt").c_str());
/*
	string cis_file(base+".cis.txt");
	ofstream out_cis(cis_file.c_str());
	string ctgs_file(base+".contigs.fasta");
	ofstream out_ctgs(ctgs_file.c_str());
*/
	bool first = true;		
	struct fa_entry* tmp = new fa_entry;
	char c = (int) in.get();
	if (c != '>') {
		cerr << "File not in fasta format\n";
		return -1;
	}
	char* scaf;
 	char buf[256];	
	in.getline(buf,256,'\n');	
	scaf = buf;
	tmp->name = "contig1";
	int width = 0;
	stringstream ss;
	int contig = 1;
	while (in.good()){
		c = (char) in.get();
		if (!in.good()) break; 
		if (c == '>'){
			print_seq(*tmp,out_ctgs);
			out_cis << tmp->name << "\t" << scaf << "\n";
			contig++;
			delete tmp;
			tmp = new fa_entry;
			in.getline(buf,256,'\n');	
			scaf = buf;
			ss << "contig" << contig;
			tmp->name = ss.str().c_str();
			ss.flush();
		} else if (is_unk(c)) {
			print_seq(*tmp,out_ctgs);
			out_cis << tmp->name << "\t" << scaf << "\n";
			contig++;
			delete tmp;
			tmp = new fa_entry;
			ss << "contig" << contig;
			tmp->name = ss.str().c_str();
			ss.flush();
			while (is_unk((char) in.peek())) in.get();	
		} else if (c == '\n') {
			continue;
		}  else {
			tmp->seq.push_front(c);
		}
	}
	print_seq(*tmp,out_ctgs);
	out_cis << tmp->name << "\t" << scaf << "\n";		
	delete tmp;
}	
