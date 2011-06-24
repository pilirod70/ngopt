#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<string>
#include<sstream>
#include<cctype>
#include<fstream>
#include<list>
#include<vector>
#include<iomanip>
#include<algorithm>
#include<set>
#include<map>
#define MAXBUFFER 100                    // max size of the gene name
#define MAX(AA,BB)  ((AA) < (BB) ? (BB) : (AA))
#define AVG(CC,DD)  (((CC)+(DD))/2)
using namespace std;

#include "point.h"			//user defined header for class point


// Global Declaration of default values and Variables



bool PRINT_BLOCK=false;
bool PRINT_B_SIMPLE=false;
bool Print_Grid= false;
bool Print_Contig= false;
bool QUIET_MODE=false;
bool TIMING=false;
string Control_File="control.txt";
string Contig_File;
string Grid_File;
string Block_File;
string Block_File_S;
bool MAX_SCORE=true;
int TOP_HITS=5;
int Max_Distance=10;
int Min_Score=200;
double Threshold=0.05;
int Min_Block_Size=3;
float BLOCK_PROB=.001;
int MaxGene=0;
unsigned int Total_Features=0;
unsigned int Total_Points=0;
float Total_Cells=0.0;
int Total_Blocks=0;
float Thresh_Dist=0;
clock_t start,finish;
int block_counter=0;  //global block counter for frequency of blocks
int MAX_BLOCK=0;        //maximum possible block size that is seen for allocating memory.
float pvalue;
float expected_blocks=0;
float p_u;   //refer documentation
float prob=0.0; /*is the probability of a cell being a point, "h" in the documentation*/
/*##############################################
#                                              #
# Declaration of data structures and typdefing #
################################################
*/


//data structure for storing gene id and direction
typedef class GENE0{
	public:
	int idx;
	int direction;
}genes;


typedef class IDENTIFYNUMBER_TO_NAME{
	public:
	int direction;
	string gene_name;
} ID2N;

vector<ID2N> id2n;     /*for retrieving the gene names and direction from identifying number*/


typedef class FEATURE_TO_GENE{
	public:
	vector<string>gene_list;
}F2G;

vector<F2G> f2g;  /* For retrieving the genes comprising each feature */


//data structure for storing gene names along with gene id and direction
//tree type data strcutre for rapid search
typedef map<string,genes> GENE_LIST;


typedef class CONTS{
	public:
	int index;
	int MAX_N;
	int feature;
} CHROMO;

vector<CHROMO> nMAX;
vector<int> nSHIFT;
vector<int> TRASH;
vector<int> FREQ;

typedef class GENE1{
		public:
		vector<int> FAMILY;
		int gene1;
		vector<CHROMO> point_id;
}FEAT;

typedef list<int> COLL;
typedef class POS{
	public:
	int gene;
	int contig;
	int feature_id;
	int flag;
	int direction;
	int where;
	COLL collapse;
}ELEM;


typedef class SUB{
	public:
	int contig1;
	int contig2;
	int points;
	float cells;
	int blocks;
	bool ext;
}SUBGRID;

vector<SUBGRID> sub_grid;

vector<FEAT> CL;

COLL::iterator datai;
vector<FEAT> LIST1;
vector<FEAT> LIST2;

typedef vector<point> ARR_POINTS;

ARR_POINTS GENE,INV_GENE, CG;

typedef class A{
	public:
	bool in_grids;
	int contig1;
	int contig2;
	int index;
	point feat;
	int avg_count;
	int flag;
	int out;
	int in;
	list<int> NH;
	int DP;
} POINTS1;

vector<POINTS1> MATCHES;
vector<POINTS1> EX;
vector<POINTS1>::iterator grid_iter;    //global iterator for type POINTS.
typedef vector<ELEM> TABLE1;

TABLE1 table(MaxGene);

typedef class B{
	public:
	list<int> block;
	bool flag;
}BLOCKS;

vector<BLOCKS> block_list;


/*data structure for storing the highest scoring nodes */
typedef class TL{
public:
int score;
int idx;
}TL1;

/*data structure used for symmetry*/

typedef map<string,int> SYMM;
typedef map<string,SYMM> MAP_LIST;
/*
#########################
#  Sorting Criterions   #
#########################
*/

/*Sort criterion for choosing the end nodes of highest scores*/
class sort_crit{
public:
bool operator()(const TL1& p1, const TL1& p2) const{

return p1.score>=p2.score;
}

};


bool sortcriterion( const point& p1, const point& p2)
{	if(p1.y<p2.y)
		return true;
	else if((p1.y==p2.y)&&(p1.x<p2.x))
		return true;
	else
		return false;
}

bool sortcriterion1( const point& p1, const point& p2)
{	if(p1.x<p2.x)
		return true;
	else if((p1.x==p2.x)&&(p1.y<p2.y))
		return true;
	else
		return false;
}

bool criterion1(const POINTS1& p1, const POINTS1& p2)
{
	if(p1.feat.y<p2.feat.y)
		return true;
	else if((p1.feat.y==p2.feat.y)&&(p1.feat.match_score>=p2.feat.match_score))
		return true;
	else
		return false;
}

bool criterion2(const POINTS1& p1, const POINTS1& p2)
{
	if(p1.feat.x<p2.feat.x)
		return true;
	else if((p1.feat.x==p2.feat.x)&&(p1.feat.match_score>=p2.feat.match_score))
		return true;
	else
		return false;
}

bool criterion(const POINTS1& p1, const POINTS1& p2)
{
	if(p1.feat.y<p2.feat.y)
		return true;
	else if((p1.feat.y==p2.feat.y)&&(p1.feat.x<p2.feat.x))
		return true;
	else
		return false;
}

bool subgrid_sort(const SUBGRID& p1, const SUBGRID& p2)
{
	if(p1.contig1<p2.contig1)
		return true;

	else if((p1.contig1==p2.contig1)&&(p1.contig2<p2.contig2))
		return true;
	else
		return false;
}

char *get_time_stamp ()
{
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  
  return asctime(timeinfo);
}


/*
##############################################################
# Declaration of other user defined header files             #
# which require the above data structures and/or variables   #
##############################################################
*/

#include "usage.h"		//header for processing command line arguments
#include "readfile.h"
#include "tandem.h"
#include "grids.h"
#include "blocks.h"
#include "dpalgo.h"
#include "output_sc.h"

/*
###############################
#    Begin of MAIN PROGRAM    #
#                             #
###############################

 */



int main(int argc, char* argv[])
{
  if(argc>1)
  process_arguments(argc, argv);
  

  double duration,dur=0,dur1;



  //Reading of Control File the source code is in the file readfile.h
  start=clock();
  read_controlfile();
  finish=clock();
  duration=double(finish-start)/CLOCKS_PER_SEC;




  //Sorting of points based on increasing y values using quick
  //sort which on an average has a complexity of n*(log n)
  start=clock();
  sort(GENE.begin(),GENE.end(),sortcriterion);
  sort(INV_GENE.begin(),INV_GENE.end(),sortcriterion1);
  finish=clock();
  dur=double(finish-start)/CLOCKS_PER_SEC;


  //Begin of detandemization the source code is in the file tandem.h
  start=clock();
  detandemize();
  finish=clock();
  dur1=double(finish-start)/CLOCKS_PER_SEC;


  //Begin of grids making , source code is in the file grids.h
  start=clock();
  grids();
  finish=clock();
  double dur3=double(finish-start)/CLOCKS_PER_SEC;




  //Begin of Block formation based on dynamic programming algorithm
  start=clock();
  blocks();
  finish=clock();
  double dur4=double(finish-start)/CLOCKS_PER_SEC;


  if(!QUIET_MODE)
  output2screen();


  if(TIMING){

  cout<<"Reading Time "<<duration<<endl;
  cout<<"Sorting Time "<<dur<<endl;
  cout<<"From genes to features "<<dur1<<endl;
  cout<<"From features to grids "<<dur3<<endl;
  cout<<"From grids to blocks "<<dur4<<endl;

  }


  return 0;
}


