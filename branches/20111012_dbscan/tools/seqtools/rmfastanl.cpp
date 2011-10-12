#include <cstdio>
#include <list>
#include <fstream>
#include <iostream>

using namespace std;

struct fa_entry {
	char* name;
	list<char> seq;
};


void print_seq(struct fa_entry& fa) {
	cout << '>' << fa.name << endl;
	list<char>::iterator it = fa.seq.begin(); 
	int len = 0;
	while (it != fa.seq.end()){
		cout << *it;
		it++;
	}
	cout << '\n';
}


int main(int argc, char** argv) {
	
	if (argc != 2) {
		cout << "Usage: rmfanl <in.fasta>\n" << "Output printed to standard out\n";
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
			print_seq(*tmp);		
			delete tmp;
			tmp = new fa_entry;
			in.getline(buf,256,'\n');	
			tmp->name=buf;
		} else if (c == '\n') {
			continue;
		}  else {
			tmp->seq.push_front(c);
		}
	}
	print_seq(*tmp);		
	delete tmp;
}	
