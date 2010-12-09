#include <cstdio>
#include <cstring>
#include <set>
#include <iostream>
#include <fstream>
using namespace std;

set<string> keepers;

void add_read(string& ln) {
	size_t pos = ln.find('\t',0);
	string tmp = ln.substr(0,pos);
	cerr << "Adding >>" << tmp << "<<" << endl;
	keepers.insert(tmp);	
}

bool keep_read(string& ln) {
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
	char line[256];


	cerr << "Found " << keepers.size() << " reads in SAM file\n";
	while (sam_in.good()){	
		sam_in.get(line,256);
		if (line[0] == '@'){ 
			cerr << "Skipping " << line << endl;
			continue;
		} else {
			string* tmp = new string(line);
			add_read(*tmp);
		}
	}
	cerr << "Found " << keepers.size() << " reads in SAM file\n";
	ifstream fq_in(argv[2]);
	string* tmp;	
	while (fq_in >> line) {
		tmp = new string(line);
		if (keep_read(*tmp)) {
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

