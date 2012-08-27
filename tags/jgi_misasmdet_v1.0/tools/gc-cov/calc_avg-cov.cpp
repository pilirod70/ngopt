#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <list>
using namespace std;

struct tally {
	double cov;
	int len;
	void operator+=(double c);
	double coverage();
};
void tally::operator+=(double c){
	cov += c;
	len++;
}

double tally::coverage(){
	return cov/((double)len);
}


int main(int argc, char** argv) {
	
	if (argc == 1) {
		cerr << "Usage: calc_avg-cov <in.plp> \n";
		return 0;
	}
	
	string plpFile = argv[1];
	ifstream plpIn(plpFile.c_str());
	list<string> order;
	map<string,tally> contig;
	string ctg;
	string str; 
	double cov;
	double avg;
	map<string,tally>::iterator it;
	int total_sites=0;
//	while (plpIn.good()){
//	while (plpIn.good()){ //;
	while(getline(plpIn,ctg,'\t')){//;
		getline(plpIn,str,'\t');
		getline(plpIn,str,'\t');
		getline(plpIn,str,'\t');
		it = contig.find(ctg);
		cov = atof(str.c_str());
		avg += cov;
		if (it == contig.end()) {
			contig[ctg].cov = cov; 
			contig[ctg].len = 1;	
			order.push_back(ctg);
//			cerr << "New contig: >" << ctg << "<"<< endl;
//			cerr << "Processed " << total_sites << " bases\n";
		}else {
			contig[ctg] += atof(str.c_str());
		}
		total_sites++;
		getline(plpIn,str);
	}
	
	cerr << "Found " << total_sites << " bases in " << contig.size() << (contig.size() >1? " contigs.\n":" contig.\n");
	fprintf(stderr,"Genome wide average coverage: %.3g\n",avg/total_sites); 
//	cerr << "Found " << total_sites << " bases in " << order.size() << " contigs.\n";
	for (list<string>::iterator ctg_it = order.begin(); ctg_it != order.end(); ctg_it++) {	
	//	cout<< *ctg_it;
		printf("%s\t%f\n",ctg_it->c_str(), contig[*ctg_it].coverage());

	}

	

}

