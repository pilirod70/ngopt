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

	while (fa_in.peek() != -1){
		c = (char) fa_in.get();
		if (c == '>'){
			cout << tmp->name << "\t" << tmp->length << endl;
			delete tmp;
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
	cout << tmp->name << "\t" << tmp->length << endl;
}	
