#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <deque>
#include <map>
using namespace std;

int main (int argc, char** argv) {
	
	istream* in;
	if (argc == 2 && strcmp(argv[1],"-h") == 0) {
		cerr << "Usage: calc_reads_gc <in.fasta/q> > <out.gc>\nread from stdin if <in.fasta/q> absent\noutput one line per read" << '\n';
		return 0;
	} else if (argc == 1) {
		in = &cin;
	} else {
		in = new ifstream(argv[1]);
	}
	char buf[1024];
	int len = 0;
	double gc = 0.0;
	bool fastq = in->peek() == '@';


	while (in->good()){
		in->getline(buf,1024);
		if (!in->good()) break;
		in->getline(buf,1024);
		len = strlen(buf);
		for (int i = 0; i < len; i++){
			switch (buf[i]) {
				case 'g':
					gc++;
					break;
				case 'G':
					gc++;
					break;
				case 'c':
					gc++;
					break;
				case 'C':
					gc++;
					break;
				case 's':
					gc++;
					break;
				case 'S':
					gc++;
					break;
				default : break;
			}
		}
		cout << gc/len << endl;
		gc = 0.0;
		if (fastq){
			in->getline(buf,1024);
			in->getline(buf,1024);
		}
	}	
	
}





