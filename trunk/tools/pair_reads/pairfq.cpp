/**
 * pairfq -- a program to re-pair fastq format sequence reads based on their illumina identifier
 * @author Andrew Tritt
 * @author Aaron Darling
 * (c) 2010, Licensed under the GPL
 * TODO: using memory-mapped I/O would probably lead to somewhat faster file processing during the writing stage
 */

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
using namespace std;

class seqloc {
public:
	seqloc() : start_offset_1(0), length_1(0), seqlength_1(0), start_offset_2(0), length_2(0), seqlength_2(0) {};
	size_t start_offset_1;
	size_t length_1;
	size_t seqlength_1;
	size_t start_offset_2;
	size_t length_2;
	size_t seqlength_2;
};

int main (int argc, char** argv) {
	
	if (argc != 2) {
		cerr << "Usage: " << argv[0] << "  <in.fastq>" << endl;
		return 0;
	}
	
	string inFile = argv[1];
	ifstream in(inFile.c_str());
	string hdr, seq, qual;

	// use an unordered_map because it is implemented with a hash function and
	// inserts and lookups are constant-time
	unordered_map<string, seqloc> pairlocs;
	string key;
	// make one pass over the file to determine the file offsets for paired reads
	// calculate some summary statistics
	size_t paircount = 0;
	size_t pairlength = 0;
	size_t unilength = 0;
	size_t start_offset = 0;
	// find out where the end of the file is
	in.seekg(-1,ios::end);
	size_t file_size = in.tellg();
	in.seekg(0,ios::beg);	// back to beginning
	size_t prev_progress = 99999;
	cout << "Reading..";
	while (start_offset < file_size){
		// print a progress message
		if( (start_offset * 100) / file_size != prev_progress ){
			prev_progress = (start_offset * 100) / file_size;
			cout << prev_progress << "%..";
			cout.flush();
		}
		// read the four components of a fastq entry
		in >> hdr;
		in >> seq;
		in >> hdr;
		in >> qual;
		// figure out where we are in the file and check for EOF
		size_t end_offset = in.tellg();	
		end_offset++; // assume 1-byte unix end-of-lines
		if(!in.good()){
			end_offset = file_size;
		}
		
		key = hdr.substr(0,hdr.length()-1);	// strip off the 1 or 2 identifier for illumina reads
		// look up the read ID and process depending on whether it's been seen before
		seqloc& sl = pairlocs[key];
		if(sl.length_1 != 0){
			sl.start_offset_2 = start_offset;
			sl.length_2 = end_offset - start_offset;
			sl.seqlength_2 = seq.length();
			paircount++;
			pairlength += sl.seqlength_1;
			pairlength += sl.seqlength_2;
			unilength -= sl.seqlength_1;
		}else{
			sl.start_offset_1 = start_offset;
			sl.length_1 = end_offset - start_offset;
			sl.seqlength_1 = seq.length();
			unilength += sl.seqlength_1;
		}
		start_offset = end_offset;
	}

	// reset for further reading below
	in.clear();

	// write summary stats
	cout << "\nPaired count: " << paircount << endl;
	cout << "Paired Mbp: " << pairlength / 1000000 << endl;
	cout << "Unpaired count: " << pairlocs.size() - paircount << endl;
	cout << "Unpaired Mbp: " << unilength / 1000000 << endl;

	// now write out read pairs
	ofstream p1out((inFile+"_p1.fastq").c_str());
	ofstream p2out((inFile+"_p2.fastq").c_str());
	ofstream upout((inFile+"_up.fastq").c_str());
	
	unordered_map<string, seqloc>::iterator it = pairlocs.begin();
	char buf[10000];	// let's hope a read never uses more than 10k chars
	cout << "Writing..";
	size_t done = 0, modulus = pairlocs.size() / 100;
	for( ; it != pairlocs.end(); it++ ){
		if(done % modulus == 0){
			cout << (100*done) / pairlocs.size() << "%..";
			cout.flush();
		}
		in.seekg(it->second.start_offset_1);	// seek to the file offset where the read lives
		in.read( buf, it->second.length_1 );	// read in the appropriate number of bytes
		buf[it->second.length_1]=0;	// NULL terminate the string
		if(it->second.length_2==0){			
			upout << buf;	// unpaired
		}else{
			p1out << buf;
			in.seekg(it->second.start_offset_2);
			in.read( buf, it->second.length_2);
			buf[it->second.length_2]=0;
			p2out << buf;
		}
		done++;
	}
	
	p1out.close();
	p2out.close();
	upout.close();
	return 0;	// signals success by unix convention
}
 


