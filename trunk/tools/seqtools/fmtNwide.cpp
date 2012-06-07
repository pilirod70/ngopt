#include <cstdio>
#include <cstdlib>
#include <list>
#include <fstream>
#include <iostream>

using namespace std;

struct fa_entry {
	char* name;
	list<char> seq;
};


void print_seq(struct fa_entry& fa, int width) {
	cout << '>' << fa.name << endl;
	list<char>::iterator it = fa.seq.begin(); 
	int len = 0;
	while (it != fa.seq.end()){
		cout << *it;
		len++;
		it++;
		if (len % width == 0)
			cout << '\n';
	}
	if (len % width != 0) 
		cout << '\n';
}


int main(int argc, char** argv) {
	
	if (argc != 3 && argc != 2) {
		cout << "Usage: fmtNwide <width> <in.fasta>\n" << "If <in.fasta> not present, read from stdin\n"<< "Output printed to standard out\n";
		return -1;
	}

	int width = atoi(argv[1]);

	istream* in; 
	if (argc == 2) {
		in = &cin;
	} else {
		in = new ifstream(argv[2]);
	}

	struct fa_entry* tmp = new fa_entry;
	char c = (int) in->get();
	if (c != '>') {
		cerr << "File not in fasta format\n";
		return -1;
	}
 	char buf[256];	
	in->getline(buf,256,'\n');	
	tmp->name=buf;
	while (in->peek() != -1){
		c = (char) in->get();
		if (c == '>'){
			print_seq(*tmp,width);		
			delete tmp;
			tmp = new fa_entry;
			in->getline(buf,256,'\n');	
			tmp->name=buf;
		} else if (c == '\n' || c == ' ') {
		}  else {
			tmp->seq.push_back(c);
		}
	}
	print_seq(*tmp,width);		
	delete tmp;
}	
