#include <cstdio>
#include <list>
#include <fstream>
#include <iostream>

using namespace std;


int main(int argc, char** argv) {
	
	if (argc != 2) {
		cout << "Usage: maxrdlen <reads_in.fasta/q>\n" << "Output printed to standard out\n";
		return -1;
	}
	bool fastq = false;
	ifstream in(argv[1]);
	char c = (char) in.get();
	if (c=='@')
		fastq = true;
	char buf[256];
	int maxlen=0;
	int len=0;
	while(in.peek() != -1){
		in.getline(buf,256);
		in.getline(buf,256);
		if (fastq){
			in.getline(buf,256);
			in.getline(buf,256);
		}
		len = strlen(buf);
		if (maxlen < len) 
			maxlen=len;
	}
	cout << maxlen << endl;
}	
