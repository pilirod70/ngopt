/*CONVERSION FROM FEATURES TO GRIDS AND TOP HITS ARE TAKEN CARE OF IN THIS FILE*/


void grids(void)
{

	FEAT temp1;
	CHROMO temp2;
	POINTS1 temp;
	int temp_var;
	int temp_score;
	vector<FEAT> CL(Total_Features);
	vector<int>::iterator pos;
	int i,j,k=0,cumulative=0,ind=0;


	if(CG.size()==0)
	{
	cerr<<"\nfish:error No Data for comparison";
	cerr<<"\nThere are no matches to be compared please check your data or refer to the\n";
	cerr<<"manual for information on the type of data which can be processed\n";
	exit(1);
	}



	for(i=0;i<CG.size();i++)
	{



		if(table[CG[i].y-1].feature_id==table[CG[i].x-1].feature_id)
			{

				if(CG[i].mark==1)
				{

					MATCHES[MATCHES.size()-1].flag=-101;

					while(sub_grid[ind].ext!=true)
					ind++;

					if(sub_grid[ind].ext==true)
					{
						sub_grid[ind].points=k-cumulative;
						cumulative=k;
					}
					ind++;

				}
				continue;
			}


		if(CL[(table[(CG[i].y)-1].feature_id)].FAMILY.size()==0)
		{
			temp.contig1=table[CG[i].y-1].contig;
			temp.contig2=table[CG[i].x-1].contig;
			temp.feat.y=table[CG[i].y-1].feature_id;
			temp.feat.x=table[CG[i].x-1].feature_id;
			temp.feat.match_score=CG[i].match_score;
			temp.feat.check=false;
			temp.feat.direction=(table[CG[i].x-1].direction)*(table[CG[i].y-1].direction);
			temp.index=k;
			temp.flag=0;
			temp.DP=0;
			temp.in=0;
			temp.in_grids=true;
			temp.avg_count=1;
			CL[table[CG[i].y-1].feature_id].gene1=k;
			CL[table[CG[i].y-1].feature_id].FAMILY.push_back(table[CG[i].x-1].feature_id);
			MATCHES.push_back(temp);
			k++;
		}
		else
		{
			pos=find(CL[table[CG[i].y-1].feature_id].FAMILY.begin(),CL[table[CG[i].y-1].feature_id].FAMILY.end(),table[CG[i].x-1].feature_id);

			if(pos==CL[table[CG[i].y-1].feature_id].FAMILY.end())
			{
				temp.contig1=table[CG[i].y-1].contig;
				temp.contig2=table[CG[i].x-1].contig;
				temp.feat.y=table[CG[i].y-1].feature_id;
				temp.feat.x=table[CG[i].x-1].feature_id;
				temp.feat.match_score=CG[i].match_score;
				temp.feat.check=false;
				temp.feat.direction=(table[CG[i].x-1].direction)*(table[CG[i].y-1].direction);
				temp.index=k;
				temp.flag=0;
				temp.DP=0;
				temp.in=0;
				temp.in_grids=true;
				temp.avg_count=1;
				CL[table[CG[i].y-1].feature_id].gene1=k;
				CL[table[CG[i].y-1].feature_id].FAMILY.push_back(table[CG[i].x-1].feature_id);
				MATCHES.push_back(temp);
				k++;

			}
			else
			{
				if(MAX_SCORE)
				{
				   temp_score=MATCHES[CL[table[CG[i].y-1].feature_id].gene1].feat.match_score;
				   MATCHES[CL[table[CG[i].y-1].feature_id].gene1].feat.match_score=(MAX(CG[i].match_score,temp_score));
				}
				else
				{
				   temp_var=MATCHES[CL[table[CG[i].y-1].feature_id].gene1].avg_count;
				   temp_score=MATCHES[CL[table[CG[i].y-1].feature_id].gene1].feat.match_score;
				   MATCHES[CL[table[CG[i].y-1].feature_id].gene1].feat.match_score=((temp_var*temp_score)+CG[i].match_score)/(temp_var+1);
				   (MATCHES[CL[table[CG[i].y-1].feature_id].gene1].avg_count)++;
				}

			}



		}


		if(CG[i].mark==1)
		{
			if(!MATCHES.empty())
			MATCHES[MATCHES.size()-1].flag=-101;

			while(sub_grid[ind].ext!=true)
			ind++;

			if(sub_grid[ind].ext==true)
			{
				sub_grid[ind].points=k-cumulative;
				cumulative=k;
			}

			ind++;
		}






	}




/*
#############################################
#                                           #
#   TOP_HITS without losing Symmetry        #
#                                           #
#############################################
*/
/* all those feature pairs are removed which are not amongst the top hits
of a given feature */

EX=MATCHES;
sort(EX.begin(),EX.end(),criterion1);

CL.clear();
CL.resize(Total_Features);
vector<FEAT> IC(Total_Features);
//vector<CHROMO>::iterator pos1;

for(i=0;i<EX.size();i++)
{
	temp2.index=EX[i].index;
	temp2.feature=EX[i].feat.x;
        CL[EX[i].feat.y].point_id.push_back(temp2);
}

sort(EX.begin(),EX.end(),criterion2);

for(i=0;i<EX.size();i++)
{
	temp2.index=EX[i].index;
	temp2.feature=EX[i].feat.y;
	IC[EX[i].feat.x].point_id.push_back(temp2);
}


for(i=0;i<CL.size();i++)
    if(!CL[i].point_id.empty())
	if(CL[i].point_id.size()>TOP_HITS)
	     for(j=TOP_HITS;j<CL[i].point_id.size();j++)
		  if(IC[CL[i].point_id[j].feature].point_id.size()>TOP_HITS)
		      for(k=TOP_HITS;k<IC[CL[i].point_id[j].feature].point_id.size();k++)
		           if(i==IC[CL[i].point_id[j].feature].point_id[k].feature)
			       {
			         TRASH.push_back(CL[i].point_id[j].index);
				 MATCHES[CL[i].point_id[j].index].in_grids=false;
			       }






EX.clear();

k=0;
ind=0;
for(grid_iter=MATCHES.begin();grid_iter!=MATCHES.end();++grid_iter)
{
   if(grid_iter->in_grids==true)
   {
     EX.push_back(*grid_iter);
     k++;
   }

   if(grid_iter->in_grids==false&&grid_iter->flag==-101&&!EX.empty())
   EX[EX.size()-1].flag=-101;

   if(!EX.empty()&&EX[EX.size()-1].flag==-101)
   {
   	while(sub_grid[ind].ext!=true)
	ind++;

	sub_grid[ind].points=k;
	k=0;
	ind++;
    }




}



MATCHES.clear();
MATCHES=EX;
EX.clear();


/* Reassigning index numbers to points after removing those points from the grid*/

for(i=0;i<MATCHES.size();i++)
MATCHES[i].index=i;


Total_Points=MATCHES.size();


/*
#############################################
#                                           #
#Counting of cells for each pair of contig  #
#############################################
*/

	for(i=0;i<sub_grid.size();i++)
	{
		if(sub_grid[i].contig1== sub_grid[i].contig2)
			sub_grid[i].cells=((nMAX[sub_grid[i].contig1-1].feature)*(nMAX[sub_grid[i].contig1-1].feature-1))/2;

		if(sub_grid[i].contig1!=sub_grid[i].contig2)
			sub_grid[i].cells=(nMAX[sub_grid[i].contig1-1].feature)*(nMAX[sub_grid[i].contig2-1].feature);


	}


/*
######################################
#                                    #
#    Printing the grids file         #
#                                    #
######################################

*/	FILE *print_grid;

	if(Print_Grid)
	{
		print_grid=fopen(Grid_File.c_str(),"w");
                fprintf(print_grid, "%s",get_time_stamp());
                fprintf(print_grid, "FISH v1.0\n");
                fprintf(print_grid, "FISH is copyright (c)2003, University of North Carolina at Chapel Hill \n");
                fprintf(print_grid, "Authors: Sugata Chakravarty and Todd J. Vision\n\n");
		fprintf(print_grid,"Control File = %s\n\n",Control_File.c_str());
		fprintf(print_grid,"-parameters\n");
		fprintf(print_grid,"\tMIN_SCORE = %d\n",Min_Score);
		fprintf(print_grid,"\tTOP_HITS = %d\n\n",TOP_HITS);
		fprintf(print_grid,"-subgrids\n");
		fprintf(print_grid,"\tcontig1\tcontig2\t  points\tcells\n");
		for(i=0;i<sub_grid.size();i++)
			fprintf(print_grid,"%10d %10d %10d %11.0f\n",sub_grid[i].contig1,sub_grid[i].contig2,sub_grid[i].points, sub_grid[i].cells);

		fprintf(print_grid,"\n\n-points\n\n");
		fprintf(print_grid,"\tpoint\tcontig1\t   contig2\tfeat1\t  feat2\t  score\n");


		for(i=0;i<MATCHES.size();i++)
			fprintf(print_grid,"%10d %10d %10d %10d %10d %10d\n",MATCHES[i].index,MATCHES[i].contig1,MATCHES[i].contig2,MATCHES[i].feat.y,MATCHES[i].feat.x, MATCHES[i].feat.match_score);

	}

	
}



