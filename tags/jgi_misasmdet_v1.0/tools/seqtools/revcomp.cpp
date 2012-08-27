#include <cstdio>
#include <cstring>
#include <list>
#include <fstream>
#include <iostream>

using namespace std;

bool fastq;

struct seq_entry {
	char* name;
	list<char> seq;
	list<char> qual;
};

char comp(char b) {
	switch(b){
		case 'A': return 'T';
		case 'a': return 't';
		case 'U': return 'A';
		case 'u': return 'a';
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

void print_seq(struct seq_entry& fa, int width) {
	if (width == -1)
		return;
	if (fastq)
		cout << '@' << fa.name << endl;
	else
		cout << '>' << fa.name << endl;
	list<char>::iterator it = fa.seq.begin(); 
	int len = 0;
	while (it != fa.seq.end()){
		cout << *it;
		len++;
		if (len % width == 0) {
			cout << '\n';
		}
		it++;
	}
	if (len % width != 0) 
		cout << '\n';
	if (fastq){
		cout << '+' << fa.name << endl;
		it = fa.qual.begin();
		len = 0;
		while (it != fa.qual.end()){
			cout << *it;
			len++;
			if (len % width == 0) {
				cout << '\n';
			}
			it++;
		}
		if (len % width != 0) 
			cout << '\n';
	}
}


int main(int argc, char** argv) {
	
//	if (argc2) {
//		cout << "Usage: revcomp <in->fasta>\n" << "Output printed to standard out\n";
//		return -1;
//	}

	istream* in;
	if (argc == 2){
		in = new ifstream(argv[1],ios::in);
	} else {
		in = &cin;
	}


	bool first = true;		
	struct seq_entry* tmp = new seq_entry;
	char c = (int) in->get();
	if (c != '>' && c != '@') {
		cerr << "File not in fasta format: " << c << "\n";
		return -1;
	}
	if (c == '@')
		fastq = true;
 	char buf[256];	
	in->getline(buf,256,'\n');	
	tmp->name = new char[strlen(buf)];
	strncpy(tmp->name,buf,strlen(buf));
	int width = 0;
	int stretch = 0;
	bool in_qual = false;
	while (in->peek() != -1){
		c = (char) in->get();
		if (c == '>' || c == '@'){
			print_seq(*tmp,width);		
			delete tmp;
			tmp = new seq_entry;
			in->getline(buf,256,'\n');	
			tmp->name = new char[strlen(buf)];
			strncpy(tmp->name,buf,strlen(buf));
			stretch = 0;
			in_qual = false;
		} else if (c == '+') {
			in->getline(buf,256,'\n');	
			in_qual = true;	
		} else if (c == '\n') {
			if (stretch > width)
				width = stretch;
			stretch = 0;
		}  else {
			if (in_qual) {
				tmp->qual.push_front(c);
			} else { 
				stretch++;
				tmp->seq.push_front(comp(c));
			}
		}
	}
	print_seq(*tmp,width);		
	delete tmp;
}	
