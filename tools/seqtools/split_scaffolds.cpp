#include <cstdio>
#include <list>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>

using namespace std;

struct fa_entry {
	char* name;
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
		cout << "Usage: rmfanl <in.fasta>\n" << "Output printed to standard out\n";
		return -1;
	}

	ifstream in(argv[1]);
	string infile(argv[1]);
	string base(infile.substr(0,infile.find_last_of(".fa")));	
	ofstream out_cis(base+".cis.txt");
	ofstream out_ctgs(base+".contigs.fasta");

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
	char itoa_buf[10];
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
			tmp->name = "contig"+buf;
		} else if (is_unk(c)) {
			print_seq(*tmp,out_ctgs);
			out_cis << tmp->name << "\t" << scaf << "\n";
			contig++;
			delete tmp;
			tmp = new fa_entry;
			tmp->name = "contig"+buf;
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
