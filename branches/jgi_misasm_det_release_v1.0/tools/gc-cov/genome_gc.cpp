#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <cstring>
using namespace std;

map<char,int> tot;

void addbase(char& b) {
	switch (b) {
		case 'a': tot['a']++; break;
		case 'A': tot['a']++; break;
		case 't': tot['t']++; break;
		case 'T': tot['t']++; break;
		case 'g': tot['g']++; break;
		case 'G': tot['g']++; break;
		case 'c': tot['c']++; break;
		case 'C': tot['c']++; break;
		case 'M': break;
		case 'm': break;
		case 'R': break;
		case 'r': break;
		case 'W': break;
		case 'w': break;
		case 'S': break;
		case 's': break;
		case 'Y': break;
		case 'y': break;
		case 'K': break;
		case 'k': break;
		case 'V': break;
		case 'v': break;
		case 'H': break;
		case 'h': break;
		case 'D': break;
		case 'd': break;
		case 'B': break;
		case 'b': break;
		case 'X': break;
		case 'x': break;
		case 'n': break;
		case 'N': break;
		default: cerr << "Unrecognizable character: " << b << endl;
	}
}

void resetCounts () {
	tot['a'] = 0;
	tot['t'] = 0;
	tot['g'] = 0;
	tot['c'] = 0;
}

double calcGC(istream& in) {

	resetCounts();
	string hdr;
	const char* split = " \t|:";
	in.get();
	getline(in,hdr);
	for(char b = in.get(); in.good(); b = in.get()){
		if (b == '>' ){ 
			getline(in,hdr);
			continue;
		} else if (b == '\n')
			continue;
		addbase(b);
	}
	
	return 100*(tot['g']+tot['c'])/(double)(tot['g']+tot['c']+tot['a']+tot['t']);
}

int main (int argc, char** argv) {
	
	if (argc == 2 && strcmp(argv[1],"-h") == 0) {
		cerr << "Usage: " << argv[0] << " <in.fasta> > <out.gc>" << '\n';
		return 0;
	}	
	double gc = 0.0;
	istream* in;
	if (argc == 1) {
		in = &cin;
		gc = calcGC(*in);
		fprintf(stdout,"%.2f\n",gc);
	} else {
		for (int i = 1; i < argc; i++) {
			in = new ifstream(argv[i]);
			gc = calcGC(*in);
			delete in;
			fprintf(stdout, "%s\t%.2f\n",argv[i],gc);
		}
	}
}

