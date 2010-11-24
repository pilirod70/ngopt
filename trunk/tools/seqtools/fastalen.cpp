#include <cstdio>
#include <list>
#include <fstream>
#include <iostream>

using namespace std;

struct fa_entry {
	char* name;
	int length;
};

int main(int argc, char** argv) {
	
	if (argc != 2) {
		cout << "Usage: fastalen <in.fa>\n" << "Output printed to standard out\n";
		return -1;
	}

	ifstream fa_in(argv[1]);

	list<struct fa_entry> seqs;

	bool first = true;		
	struct fa_entry* tmp = new fa_entry;
	char c = (int) fa_in.get();
	if (c != '>') {
		cerr << "File not in fasta format\n";
		return -1;
	}
 	char buf[256];	
	fa_in.getline(buf,256,'\n');	
	tmp->name=buf;
	tmp->length = 0;	

	while (fa_in.good()){
		c = (char) fa_in.get();
		if (c == '>'){
			seqs.push_back(*tmp);
			tmp = new fa_entry;
			fa_in.getline(buf,256,'\n');	
			tmp->name=buf;
			tmp->length = 0;			
		} else if (c == '\n') {
			continue;
		}  else {
			tmp->length++;
		}
	}
	seqs.push_back(*tmp);
	for(list<struct fa_entry>::iterator it = seqs.begin();it!=seqs.end();it++){
		cout << it->name << "\t" << it->length << endl;
	}
}	
