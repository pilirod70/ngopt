#include <cstdio>
#include <list>
#include <fstream>
#include <iostream>

using namespace std;

struct fa_entry {
	char* name;
	list<char> seq;
};

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

void print_seq(struct fa_entry& fa, int width) {
	cout << '>' << fa.name << endl;
	list<char>::iterator it = fa.seq.begin(); 
	int len = 0;
	while (it != fa.seq.end()){
		cout << *it;
		len++;
		it++;
		if (len % width == 0) {
			cout << '\n';
		}
	}
		
	if (len % width != 0) {
		cout << '\n';
	}
}


int main(int argc, char** argv) {
	
	if (argc != 2) {
		cout << "Usage: revcomp <in.fasta>\n" << "Output printed to standard out\n";
		return -1;
	}

	ifstream in(argv[1]);


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
			print_seq(*tmp,width);		
			delete tmp;
			tmp = new fa_entry;
			in.getline(buf,256,'\n');	
			tmp->name=buf;
			stretch = 0;
		} else if (c == '\n') {
			if (stretch > width)
				width = stretch;
			stretch = 0;
		}  else {
			stretch++;
			tmp->seq.push_front(comp(c));
		}
	}
	print_seq(*tmp,width);		
	delete tmp;
}	
