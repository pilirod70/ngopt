#include <cstdio>
#include <cstring>
#include <list>
#include <fstream>
#include <iostream>

using namespace std;

struct fq_entry {
	char* name;
	list<char> seq;
	list<char> qual;
};


void print_seq(struct fq_entry& fa) {
	cout << '@' << fa.name << endl;
	list<char>::iterator it = fa.seq.begin(); 
	int len = 0;
	while (it != fa.seq.end()){
		cout << *it;
		it++;
	}
	cout << '\n';
	cout << '+' << fa.name << endl;
	it = fa.qual.begin();
	while (it != fa.qual.end()){
		cout << *it;
		it++;
	}
	cout << '\n';
}


int main(int argc, char** argv) {
	
	if (argc != 2) {
		cout << "Usage: rmfqnl <in.fastq>\n" << "Output printed to standard out\n";
		return -1;
	}

	ifstream in(argv[1]);


	bool first = true;		
	struct fq_entry* tmp = new fq_entry;
	char c = (int) in.get();
	if (c != '@') {
		cerr << "File not in fastq format\n";
		return -1;
	}
 	char buf[256];	
	in.getline(buf,256,'\n');	
	tmp->name = new char[strlen(buf)];
	strcpy(tmp->name,buf);
	bool in_qual = false;
	while (in.peek() != -1){
		c = (char) in.get();
		if (c == '@'){
			print_seq(*tmp);		
			in_qual = false;
			delete tmp;
			tmp = new fq_entry;
			in.getline(buf,256,'\n');	
			tmp->name = new char[strlen(buf)];
			strcpy(tmp->name,buf);
		} else if (c == '\n') {
			continue;
		} else if (c == '+') {
			in.getline(buf,256,'\n');	
			in_qual = true;	
		}  else {
			if (in_qual) {
				tmp->qual.push_front(c);
			} else { 
				tmp->seq.push_front(c);
			}
		}
	}
	print_seq(*tmp);		
	delete tmp;
}	
