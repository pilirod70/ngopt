/*##################################
## Tandemization of genes	  ##
##				  ##
##				  ##
####################################
*/


void detandemize(void)
{
	int i;
	vector<int>::iterator pos;

	FEAT temp1;	//temporary variable of type FEAT for declaration of FEAT see "fish.cpp"

	/*The following "for" loop forms a array of lists, each row contains all the matches of a
	particular gene. row1: 1-->(3,4,5,6)  row2: 3-->(4,6,7) etc.*/

	for( i=0;i<GENE.size();i++)
	{

		temp1.gene1=GENE[i].y;

		if (temp1.FAMILY.size()==0)
			temp1.FAMILY.push_back(GENE[i].y);

		temp1.FAMILY.push_back(GENE[i].x);

		if((i==GENE.size()-1)||(GENE[i+1].y>GENE[i].y))
		{
			pos=temp1.FAMILY.begin();
			pos++;
			sort(pos,temp1.FAMILY.end());
			LIST1.push_back(temp1);
			temp1.FAMILY.clear();

		}



	}


	temp1.FAMILY.clear();

	for( i=0;i<INV_GENE.size();i++)
	{

		temp1.gene1=INV_GENE[i].x;

		if (temp1.FAMILY.size()==0)
			temp1.FAMILY.push_back(INV_GENE[i].x);

		temp1.FAMILY.push_back(INV_GENE[i].y);

		if((i==INV_GENE.size()-1)||(INV_GENE[i+1].x>INV_GENE[i].x))
		{
			pos=temp1.FAMILY.begin();
			pos++;
			sort(pos,temp1.FAMILY.end());
			LIST2.push_back(temp1);
			temp1.FAMILY.clear();

		}



	}


 table.resize(MaxGene);


int p=0,count=0;

/*intitalizing the entries of the vector "table" which contains entries of type "ELEM" */
for( i=0;i<MaxGene;i++)
{
	table[i].flag=0;
	table[i].where=0;
	table[i].gene=i+1;
	table[i].feature_id=-10;

		if(i+1-nSHIFT[count]<nMAX[count].MAX_N)
			table[i].contig=nMAX[count].index;



		if(i+1-nSHIFT[count]==nMAX[count].MAX_N)
		{
			table[i].contig=nMAX[count].index;
			count++;
		}



}



int c_index=0;

for(i=0;i<LIST1.size();i++)
{
	table[LIST1[i].gene1-1].flag=1;
	int j=0;
	c_index=LIST1[i].FAMILY[j];
	while(j+1<LIST1[i].FAMILY.size())
	{

		if(LIST1[i].FAMILY[j+1]-c_index<=Max_Distance)
		{

			table[LIST1[i].FAMILY[j+1]-1].where=c_index;
			table[c_index-1].collapse.insert(table[c_index-1].collapse.end(),LIST1[i].FAMILY[j+1]);
		}
		else
			c_index=LIST1[i].FAMILY[j+1];


		j++;
	}
}

c_index=0;


for(i=0;i<LIST2.size();i++)
{
	table[LIST2[i].gene1-1].flag=1;
	int j=1;
	c_index=LIST2[i].FAMILY[j];
	while(j+1<LIST2[i].FAMILY.size())
	{

		if(LIST2[i].FAMILY[j+1]-c_index<=Max_Distance)
		{

			table[LIST2[i].FAMILY[j+1]-1].where=c_index;
			table[c_index-1].collapse.insert(table[c_index-1].collapse.end(),LIST2[i].FAMILY[j+1]);
		}
		else
			c_index=LIST2[i].FAMILY[j+1];


		j++;
	}
}


	cout<<endl;
	for(i=table.size()-1;i>=0;i--)
	{
		if(table[i].where!=0)
		{
			//cout<<table[i].gene<<"    "<<table[i].where<<endl;
			table[table[i].where-1].collapse.insert(table[table[i].where-1].collapse.end(),table[i].collapse.begin(),table[i].collapse.end());
			table[i].collapse.clear();
			table[i].flag=-1;
		}
	}

int temp_dir=0;

/*loop for counting the total number of features for each contig */
/*Loop for assigning directions to features */

	for(i=0;i<table.size();i++)
	{
		if(table[i].collapse.size()>0)
		{
			table[i].collapse.sort();
			table[i].collapse.unique();
		}


		if(table[i].flag==0||table[i].collapse.size()>0||table[i].flag==1)
			nMAX[table[i].contig-1].feature=nMAX[table[i].contig-1].feature+1;




		if(table[i].flag==0)
		table[i].direction=id2n[i].direction;

		if(table[i].collapse.size()>0)
		{
		  temp_dir=id2n[i].direction;

		  for(datai=table[i].collapse.begin();datai!=table[i].collapse.end();++datai)
		  temp_dir=temp_dir+ id2n[*datai-1].direction;

		  if(temp_dir>0)
		  table[i].direction=1;

		  else if(temp_dir<0)
		  table[i].direction=-1;

		  else if(temp_dir==0)
		  table[i].direction=0;

		}

		if(table[i].flag==1 && table[i].where==0)
		table[i].direction=id2n[i].direction;

		if(table[i].flag==-1)
		table[i].direction=table[table[i].where-1].direction;


	}







/*
################################
#Begin Printing the contig file#
#and feature_id allocation     #
################################
*/

	ofstream contig_print;

	if(Print_Contig)
	{

		contig_print.open(Contig_File.c_str(),ios::out);

		contig_print<<get_time_stamp()<<"";
                contig_print<<"FISH v1.0\n";
                contig_print<<"FISH is copyright (c)2003, University of North Carolina at Chapel Hill\n";
                contig_print<<"Authors: Sugata Chakravarty and Todd J. Vision\n\n";
		contig_print<<"-parameters\n\tMIN_SCORE = "<<Min_Score<<"\n\tMAX_DIST = "<<Max_Distance<<"\n\n";
		contig_print<<"-contigs\n\n";
		contig_print<<setw(10)<<"contig"<<setw(10)<<"markers"<<setw(10)<<"features"<<endl;

		for(i=0;i<nMAX.size();i++)
			contig_print<<setw(10)<<nMAX[i].index<<setw(10)<<nMAX[i].MAX_N<<setw(10)<<nMAX[i].feature<<endl;

		contig_print<<"\n\n\n";
		contig_print<<"-features\n\n";
		contig_print<<setw(10)<<"feature"<<setw(10)<<"contig"<<setw(12)<<"orientation"<<setw(10)<<"genes"<<endl;
	}

		int k=0,index=0,prev=0;
		for(i=0;i<table.size();i++)
		{

			if(table[i].flag==0&&table[i].collapse.size()==0)
			{
				if(Print_Contig)
				contig_print<<setw(10)<<k<<setw(10)<<table[i].contig<<setw(12)<<table[i].direction<<setw(10)<<id2n[i].gene_name<<"\n";
				table[i].feature_id=k;
				k++;
				continue;
			}

			if(table[i].flag==1&&table[i].collapse.size()==0)
			{
				if(Print_Contig)
				contig_print<<setw(10)<<k<<setw(10)<<table[i].contig<<setw(12)<<table[i].direction<<setw(10)<<id2n[i].gene_name<<"\n";
				table[i].feature_id=k;
				k++;
				continue;
			}

			if(table[i].collapse.size()>0)

			{
				if(Print_Contig)
				contig_print<<setw(10)<<k<<setw(10)<<table[i].contig<<setw(12)<<table[i].direction<<setw(10)<<id2n[i].gene_name<<"\n";
				table[i].feature_id=k;


				for(datai=table[i].collapse.begin();datai!=table[i].collapse.end();++datai)
				{
					table[*datai-1].feature_id=k;
					if(Print_Contig)
					contig_print<<setw(10)<<k<<setw(10)<<table[i].contig<<setw(12)<<table[i].direction<<setw(10)<<id2n[*datai-1].gene_name<<"\n";

				}
				k++;
				continue;
			}



		}

		
		Total_Features=k;

		f2g.resize(Total_Features);

		for(i=0;i<table.size();i++)
		f2g[table[i].feature_id].gene_list.push_back(id2n[i].gene_name);

/*################################
## Free the memory By deleting  ##
## GENE data structure and LIST ##
##################################
*/
		//INV_GENE.clear();
		//GENE.clear();
		LIST1.clear();
		LIST2.clear();





}








