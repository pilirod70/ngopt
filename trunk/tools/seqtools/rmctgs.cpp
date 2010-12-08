#include <cstdio>
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
	
	if (argc != 3) {
		cout << "Usage: rmctgs <minlen> <in.fasta>\n" << "Output printed to standard out\n";
		return -1;
	}

	int minlen = atoi(argv[1]);
	ifstream in(argv[2]);


	bool first = true;		
	struct fa_entry* tmp = new fa_entry;
	char c = (int) in.get();
	if (c != '>') {
		cerr << "File not in fasta format\n";
		return -1;
	}
 	char buf[256];	
	in.getline(buf,256,'\n');	
	tmp->name=buf;
	int width = 0;
	int stretch = 0;
	while (in.peek() != -1){
		c = (char) in.get();
		if (c == '>'){
			if (tmp->seq.size() >= minlen) 
				print_seq(*tmp,width);		
			delete tmp;
			tmp = new fa_entry;
			in.getline(buf,256,'\n');	
			tmp->name=buf;
			stretch = 0;
			width = 0;
		} else if (c == '\n') {
			if (stretch > width)
				width = stretch;
			stretch = 0;
		}  else {
			stretch++;
			tmp->seq.push_back(c);
		}
	}
	if (tmp->seq.size() >= minlen) 
		print_seq(*tmp,width);		
	delete tmp;
}	
