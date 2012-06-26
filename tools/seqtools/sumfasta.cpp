#include <cstdio>
#include <list>
#include <fstream>
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

using namespace std;

void printStats(char* file) {
	ifstream fa_in(file);

	bool first = true;		
	char c = (char) fa_in.get();
	if (c != '>') {
		cerr << "File " << file << " not in fasta format\n" << "c = " << c << endl;
		return;
	}
	int tot = 0;
	int numN = 0;
	int numNuc = 0;
	int len = 0;
 	char buf[256];	
	fa_in.getline(buf,256,'\n');	
	vector<int> data;
	while (fa_in.peek() != -1){
		c = (char) fa_in.get();
		if (c == '>'){
			data.push_back(len);
			len = 0;
			fa_in.getline(buf,256,'\n');
		} else if (c == '\n') {
			continue;
		}  else {
			if (c == 'n' || c == 'N') {
				numN++;
			} else {
				numNuc++;
			}
			len++;
		}
	}
	data.push_back(len);
	tot = numN + numNuc;
	sort (data.begin(),data.end());
	int tally = 0;
	int n50 = 0;
	int half = tot/2;
	for (vector<int>::iterator it = data.begin(); it != data.end(); it++){
		tally += *it;
		if (tally > half) {
			n50 = *it;
			break;
		}
	}
	// Fasta_File	Num_Ctgs	N50	AvgCtgLen	MaxCtgLen	MinCtgLen	Total_Bases	Num_Nucs	Num_Unks
	cout << file << "\t" << data.size() << "\t" << n50 << "\t" << (tot/data.size()) << "\t" << data.back() << "\t" << *data.begin() << "\t" << tot << "\t" << numNuc << "\t" << numN << endl;
}

int main(int argc, char** argv) {
	
	if (argc != 2) {
		cout << "Usage: sumfasta <fasta1> <fasta2> ... <fastaN>\n" << "Output printed to standard out\n";
		return -1;
	}
	for (int i = 1; i < argc; i++){
		printStats(argv[i]);
	}

}	
