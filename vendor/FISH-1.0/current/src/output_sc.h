void output2screen(void)
{  
   int i,j;

   printf("%s",get_time_stamp());
   printf("FISH v1.0\n");
   printf("FISH is copyright ©2003, University of North Carolina at Chapel Hill \n");
   printf("Authors: Sugata Chakravarty and Todd J. Vision\n\n");
   printf("Control File = %s\n\n",Control_File.c_str());
   printf("processing contigs\n");
   printf("\tmin_score = %d\n\tmax_dist = %d\n\n",Min_Score,Max_Distance);

   cout<<setw(10)<<"contig"<<setw(10)<<"msrkers"<<setw(10)<<"features"<<endl;

   for(i=0;i<nMAX.size();i++)
   cout<<setw(10)<<nMAX[i].index<<setw(10)<<nMAX[i].MAX_N<<setw(10)<<nMAX[i].feature<<endl;

   if(Print_Contig)
   printf("\nwriting contig output to %s \n\n",Contig_File.c_str());

   printf("\nprocessing matches\n");
   printf("\tmin_score = %d\n\ttop_hits = %d\n\n",Min_Score,TOP_HITS);
   printf("\tcontig1\tcontig2\t  points       cells\n\n");
   for(i=0;i<sub_grid.size();i++)
   printf("%10d %10d %10d %12.0f\n",sub_grid[i].contig1,sub_grid[i].contig2,sub_grid[i].points, sub_grid[i].cells);

   if(Print_Grid)
   printf("\nwriting grid output to %s\n\n",Grid_File.c_str());

   printf("\nprocessing blocks\n");
   printf("\tT = %0.2f\n",Threshold);
   printf("\ttotal_points = %d\n",Total_Points);
   printf("\ttotal_features = %d\n",Total_Features);
   printf("\ttotal_cells = %9.0f\n",Total_Cells);
   printf("\tprobability h = %f\n",prob);
   printf("\td_T = %f\n",Thresh_Dist);
   printf("\tmin_edges = %d\n",Min_Block_Size-1);
   printf("\tblock_prob = %1.5f\n",BLOCK_PROB);
   printf("\ttotal_blocks = %d\n\n",Total_Blocks);

   cout<<setw(10)<<"contig1"<<setw(10)<<"contig2"<<setw(10)<<"blocks\n";
   for(i=0;i<sub_grid.size();i++)
   cout<<setw(9)<<sub_grid[i].contig1<<setw(9)<<sub_grid[i].contig2<<setw(9)<<sub_grid[i].blocks<<endl;

   if(PRINT_BLOCK)
   printf("\nwriting output to %s\n",Block_File.c_str());

   printf("\nblock statistics\n\n");
   cout<<"points\tp-value\t   frequency"<<endl;


   int d_t,n;
   d_t=int(Thresh_Dist);
   n=d_t*(d_t-1);
   float temp_prob;


   for(i=0;i<FREQ.size();i++)
   {
     if(FREQ[i]>0)
        {
		p_u=pow(n*prob,i)*prob;
		pvalue=1-exp(-1*Total_Cells*p_u);
		expected_blocks=Total_Cells*p_u;
		printf("%d\t%.2e\t|",i+1,pvalue);

		if(FREQ[i]<=50)
		{
		temp_prob=expected_blocks/FREQ[i];

		for(j=0;j<FREQ[i];j++)
		printf("*");

		printf("%d",FREQ[i]);

		if(temp_prob>BLOCK_PROB)
		printf("!\n");
		else
		printf("\n");


		}
		else
		{
			temp_prob=expected_blocks/FREQ[i];

			for(j=0;j<47;j++)
			printf("*");

			printf("...%d",FREQ[i]);

			if(temp_prob>BLOCK_PROB)
			printf("!\n");
		else
			printf("\n");


		}

         }
    }




     printf(" \t \t \t|----|----|----|----|----|----|----|----|----|----|\n");
     printf(" \t \t \t0         10        20        30        40        50\n");












}


