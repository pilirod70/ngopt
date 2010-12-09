#include <cstdio>
#include <cstring>
#include <unordered_set>
#include <iostream>
#include <fstream>
using namespace std;


unordered_set<string> keepers;

void add_read(string ln) {
	size_t pos = ln.find('\t',0);
	keepers.insert(ln.substr(0,pos));	
}

bool keep_read(string ln) {
	size_t pos = ln.find('/',0);
	return keepers.find(ln.substr(1,ln.length()-2)) != keepers.end();
}


void usage() {
	cout << "Usage: exmap <mapped.sam> <reads.fastq>\n";
	cout << "Output is printed to standard output.\n";
}

int main (int argc, char** argv) {
	
	if (argc != 3) {
		usage();	
	}

	ifstream sam_in(argv[1]);
		
	string line;
	while (sam_in >> line){
		if (line[0] == '@') 
			continue;
		else 
			add_read(line);
	}
	
	ifstream fq_in(argv[2]);
	
	while (fq_in >> line) {
		if (keep_read(line)) {
			cout << line << endl;
			fq_in >> line;
			cout << line << endl;
			fq_in >> line;
			cout << line << endl;
			fq_in >> line;
			cout << line << endl;
		} else {
			fq_in >> line;
			fq_in >> line;
			fq_in >> line;
		}
	}
	
	

}

