void dpalgo(void);
void blocks(void)
{
	/*Calculation of Threshold distance */

	int i,j,l;
	int d_t,n;
	block_counter=0;
	vector<POINTS1>:: iterator pos;
	vector<POINTS1>:: iterator pos1;
	list<int>:: iterator dataj;
	for(i=0;i<sub_grid.size();i++)
		Total_Cells=Total_Cells+sub_grid[i].cells;


	prob= float(Total_Points)/Total_Cells;
	Thresh_Dist =sqrt((float((log(1-Threshold)))/(log(1-prob)))+.25) + .5;
	d_t=int(Thresh_Dist);
	n=d_t*(d_t-1);


	if(MATCHES.size()<2)
	{
		cerr<<"No matches to be compared. Please check your data\n";
		exit(1);
	}

	pos1=MATCHES.begin();
	pos=MATCHES.begin();
	for(i=0;i<sub_grid.size();i++)
	{
		if(sub_grid[i].points==0)
		continue;

		pos1=pos1+sub_grid[i].points-1;
		EX.assign(pos,pos1+1);
		sort(EX.begin(),EX.end(),criterion);
		pos=++pos1;


		dpalgo();
		sub_grid[i].blocks=block_counter;
		block_counter=0;

		if(block_list.size() >0)
                {
		  block_list[block_list.size()-1].flag=true;
                }

		EX.clear();
	}



	FREQ.resize(MAX_BLOCK);
	for(i=0;i<FREQ.size();i++)
	FREQ[i]=0;


	for(i=0;i<block_list.size();i++)
	FREQ[block_list[i].block.size()-1]=FREQ[block_list[i].block.size()-1]+1;

Total_Blocks=block_list.size();
if(PRINT_BLOCK)
{


    int tempx, tempy;
    int tempdist=0;
    ofstream p2b;
    p2b.open(Block_File.c_str(),ios::out);
    p2b<<get_time_stamp()<<"";
    p2b<<"FISH v1.0\n";
    p2b<<"FISH is copyright (c)2003, University of North Carolina at Chapel Hill\n";
    p2b<<"Authors: Sugata Chakravarty and Todd J. Vision\n\n";
    p2b<<"Control File = "<<Control_File<<"\n\n";
    p2b<<"-parameters\n";
    p2b<<"\tT = "<<Threshold<<"\n\ttotal points = "<<Total_Points<<"\n";
    p2b<<"\ttotal_features = "<<Total_Features<<"\n\ttotal_cells = "<<Total_Cells<<"\n";
    p2b<<"\th = "<<prob<<"\n\td_T = "<<Thresh_Dist<<"\n\tMin_Block_Size = "<<Min_Block_Size;
    p2b<<"\n\tblock_prob = "<<BLOCK_PROB<<"\n\ttotal_blocks = "<<block_list.size()<<"\n\n\n";

    for(i=0;i<block_list.size();i++)
    {

      dataj=block_list[i].block.begin();
      p2b<<"-block "<<i<<endl;
      p2b<<setw(10)<<"points"<<setw(10)<<"contig1"<<setw(10)<<"contig2"<<endl;
      p2b<<setw(9)<<block_list[i].block.size()<<setw(9)<<MATCHES[*dataj].contig1<<setw(9)<<MATCHES[*dataj].contig2<<"\n\n";
      p2b<<"\tpoint\tdist\tmarkers\torientation"<<endl;

      tempy=MATCHES[*dataj].feat.y;
      tempx=MATCHES[*dataj].feat.x;

      for(dataj=block_list[i].block.begin();dataj!=block_list[i].block.end();++dataj)
	{
	  tempdist=(MATCHES[*dataj].feat.y-tempy)+abs(MATCHES[*dataj].feat.x-tempx);

	  p2b<<"\t"<<*dataj<<"\t"<<tempdist<<"\t{";

	  for(j=0;j<f2g[MATCHES[*dataj].feat.y].gene_list.size();j++)
	  p2b<<f2g[MATCHES[*dataj].feat.y].gene_list[j]<<" ";

	  p2b<<"}{";

	  for(j=0;j<f2g[MATCHES[*dataj].feat.x].gene_list.size();j++)
	  p2b<<f2g[MATCHES[*dataj].feat.x].gene_list[j]<<" ";

	  p2b<<"}\t"<<MATCHES[*dataj].feat.direction<<"\n";
	  tempx=MATCHES[*dataj].feat.x;
	  tempy=MATCHES[*dataj].feat.y;

	}


      p2b<<"\n\n";

    }


    p2b<<"-by contig\n";
    p2b<<setw(10)<<"contig1"<<setw(10)<<"contig2"<<setw(10)<<"blocks\n";
    for(i=0;i<sub_grid.size();i++)
    p2b<<setw(9)<<sub_grid[i].contig1<<setw(9)<<sub_grid[i].contig2<<setw(9)<<sub_grid[i].blocks<<endl;

    p2b<<"\n-by size\n";
    p2b<<setw(10)<<"points"<<setw(10)<<"obs"<<setw(10)<<"exp"<<setw(25)<<"p\n";
    for(l=0;l<FREQ.size();l++)
    {
      if(FREQ[l]>0)
        {
	  p_u=pow(n*prob,l)*prob;
          pvalue=1-exp(-1*Total_Cells*p_u);
	  expected_blocks=Total_Cells*p_u;
          p2b<<setw(9)<<l+1<<setw(9)<<FREQ[l]<<setw(15)<<scientific<<setprecision(2)<<expected_blocks<<setw(25)<<scientific<<pvalue<<"\n";
	}
    }


}




if(PRINT_B_SIMPLE)
{
	ofstream block_print;

	block_print.open(Block_File_S.c_str(),ios::out);


        
	block_print<<block_list.size()<<endl;
	block_print<<MATCHES.size()<<"\n";

	for(i=0;i<block_list.size();i++)
	{
		block_print<<block_list[i].block.size();

		//if(block_list[i].block.size()>=Min_Block_Size)
		{
			for(dataj=block_list[i].block.begin();dataj!=block_list[i].block.end();++dataj)
			{
				block_print<<","<<*dataj;
			}

			block_print<<"*"<<endl;
		}
	}

	block_print.close();
}



}




