#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <map>
#include <list>
using namespace std;

struct site {
	string ctg;
	int pos;
	char b;
	int cov;
	double gc;
	bool complete(){ return cov > 0 && gc >= 0.0; };
};



void printSite(site& s, ostream& out){
	out << s.ctg <<'\t'<< s.pos <<'\t'<< s.b <<'\t'<< s.cov <<'\t'<< s.gc << endl;
}

list<string> order; 
map<string,map<int,site> > table;
bool KEEP_ALL = false;

int main (int argc, char** argv) {
	
//	cerr << "Before" << endl;
	if (argc == 1){
		cerr << "Usage: merge_gc-plp " << "<in.plp> <in.gc> > <out.cvgc>\n";
		return 0;  
	}
//	cerr << "After" << endl;

	string str = argv[1];
	ifstream cov_in(str.c_str());
	string ctg;
	int pos;
	char b;
	int cov;
	site* s;
	map<int,site>* sites;
	while (cov_in.good()){
		getline(cov_in,ctg,'\t');
		getline(cov_in,str,'\t');
		pos = atoi(str.c_str());
		getline(cov_in,str,'\t');
		b = str.at(0);
		getline(cov_in,str,'\t');
		cov = atoi(str.c_str());
		getline(cov_in,str);
		s = new site; 
		s->ctg = ctg;
		s->pos = pos;
		s->b = b;
		s->cov = cov;
		s->gc = -1.0;
		sites = &(table.find(ctg)->second);
		if (sites == &(table.end()->second)){
			sites = new map<int,site>();
			(*sites)[pos] = *s;	
			table[ctg] = *sites;
			order.push_back(ctg);
		} else {
			(*sites)[pos] = *s;	
		}
	}	
	str = argv[2];
	ifstream gc_in(str.c_str());
	double gc;
//	while(gc_in.good()){
	while(getline(gc_in,ctg,'\t')){ //;
		getline(gc_in,str,'\t');
		pos = atoi(str.c_str());
		getline(gc_in,str);
		gc = atof(str.c_str());
		sites = &(table.find(ctg)->second);
		if (sites == &(table.end()->second)){
			sites = new map<int,site>();
			s = new site; 
			s->ctg = ctg;
			s->pos = pos;
			s->b = '-';
			s->cov = 0;
			s->gc = gc;
			(*sites)[pos] = *s;
			table[ctg] = *sites;
			order.push_back(ctg);
		} else {
			s = &(sites->find(pos)->second);
			if (s == &(sites->end()->second)){
				s = new site; 
				s->ctg = ctg;
				s->pos = pos;
				s->b = '-';
				s->cov = 0;
				s->gc = gc;
				(*sites)[pos] = *s;
			} else 
				s->gc = gc;
		}
	}

	list<string>::iterator ctg_it;
	map<int,site>::iterator site_it;
	cout << "ctg\tpos\tbase\tcov\tgc" << endl; 
	for (ctg_it = order.begin(); ctg_it != order.end(); ctg_it++){
		sites = &(table.find(*ctg_it)->second);
		for (site_it = sites->begin(); site_it != sites->end(); site_it++){
			s = &(site_it->second);
			if (s->complete() || KEEP_ALL)
				printSite(*s,cout);
		}	
	}
	

}








