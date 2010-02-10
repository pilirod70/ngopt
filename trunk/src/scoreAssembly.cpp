#include <iostream>
#include <libMems/IntervalList.h>
#include <libMems/MatchList.h>
#include "libMems/Backbone.h"
using namespace std;
using namespace mems;
using namespace genome;

void tallyGaps(vector<string>& aln, vector<int>& gaps, int seq);

int main( int argc, char* argv[] )
{
	if( argc != 2 )
	{
		cerr << "Usage: scoreAssembly <input xmfa>\n";
		return -1;
	}
	ifstream aln_in;
	aln_in.open( argv[1] );
	if( !aln_in.is_open() ){
		cerr << "Error opening " << argv[1] << endl;
		return -1;
	}
	IntervalList iv_list;
	iv_list.ReadStandardAlignment( aln_in );
	aln_in.close();
	cout << "Read alignment with " << iv_list.size() << " intervals\n";
	LoadSequences(iv_list, &cout);

	// count excess contigs
	cout << "reference replicons:\t" << iv_list.seq_table[0]->contigListSize() << endl;
	cout << "assembly contigs:\t" << iv_list.seq_table[1]->contigListSize() << endl;

	// count breakpoints in the alignment
	vector< gnSeqI > bps;
	vector<AbstractMatch*> iv_matches;
	for(int ivI=0; ivI<iv_list.size(); ivI++){
		iv_matches.push_back(&iv_list[ivI]);
	}
	IdentifyBreakpoints(iv_matches, bps);
	cout << "breakpoints\t" << bps.size() << endl;

	// count substitutions in the alignment
	int submap[128][128];
	int subtotal=0;
	vector<string> aln(2);
	for(int ivI=0; ivI < iv_list.size(); ivI++){
		GetAlignment(iv_list[ivI], iv_list.seq_table, aln);
		for(int i=0; i<aln[0].size(); i++){
			if(aln[0][i]=='-'||aln[1][i]=='-')	continue;
			submap[toupper(aln[1][i])][toupper(aln[0][i])]++;
			if(toupper(aln[0][i])!=toupper(aln[1][i])){
				subtotal++;
			}
		}
	}
	cout << "A->A\t" << submap['A']['A'] << endl;
	cout << "A->C\t" << submap['C']['A'] << endl;
	cout << "A->G\t" << submap['G']['A'] << endl;
	cout << "A->T\t" << submap['T']['A'] << endl;

	cout << "C->A\t" << submap['A']['C'] << endl;
	cout << "C->C\t" << submap['C']['C'] << endl;
	cout << "C->G\t" << submap['G']['C'] << endl;
	cout << "C->T\t" << submap['T']['C'] << endl;

	cout << "G->A\t" << submap['A']['G'] << endl;
	cout << "G->C\t" << submap['C']['G'] << endl;
	cout << "G->G\t" << submap['G']['G'] << endl;
	cout << "G->T\t" << submap['T']['G'] << endl;

	cout << "T->A\t" << submap['A']['T'] << endl;
	cout << "T->C\t" << submap['C']['T'] << endl;
	cout << "T->G\t" << submap['G']['T'] << endl;
	cout << "T->T\t" << submap['T']['T'] << endl;

	// count gaps in the alignment
	vector<int> ref_gaps;
	vector<int> ass_gaps;
	for(int ivI=0; ivI < iv_list.size(); ivI++){
		GetAlignment(iv_list[ivI], iv_list.seq_table, aln);
		tallyGaps(aln, ref_gaps, 0);
		tallyGaps(aln, ass_gaps, 1);
	}
	cout << "Gaps in reference:\n";
	for(int i=0; i<ref_gaps.size(); i++)
		cout << ref_gaps[i] << endl;
	cout << "Gaps in assembly:\n";
	for(int i=0; i<ass_gaps.size(); i++)
		cout << ass_gaps[i] << endl;
	return 0;
}

void tallyGaps(vector<string>& aln, vector<int>& gaps, int seq){
	bool gap=aln[seq][0]=='-';
	size_t gstart = 0;
	for(int i=1; i<aln[seq].size(); i++){
		if(gap&&aln[seq][i]!='-')
		{
			gaps.push_back(i-gstart);
			gap=false;
		}else if(!gap&&aln[seq][i]=='-'){
			gstart=i;
			gap=true;
		}
	}
}
