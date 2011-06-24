void usage(void);
void check(int,char* [], int);
void process_arguments(int argc, char* argv[])
{
  int i=1;string temp=argv[1];

  while(i<=argc-1)
   {
     temp=argv[i];
	  if(temp=="-f")
	    {
	    check(i,argv,argc);
	    Control_File=argv[++i];

	    }

	  else if(temp=="-h")
              usage();


	  else if(temp=="-b")
	    {
	      check(i,argv,argc);
	      Block_File=argv[++i];
	      PRINT_BLOCK=true;

	    }
	    
	    else if(temp=="-B")
	    {
	      check(i,argv,argc);
	      Block_File_S=argv[++i];
	      PRINT_B_SIMPLE=true;

	    }

	  else if(temp=="-C")
	   {
	      check(i,argv,argc);
	      Contig_File=argv[++i];
	      Print_Contig=true;

	   }

	  else if(temp=="-off")
	    Max_Distance=0;
	    
         else if(temp=="-q")
	    QUIET_MODE=true;
	    
	 else if(temp=="-t")
	 TIMING=true;


	  else if(temp=="-g")
	    {
	      check(i,argv,argc);
	      Grid_File=argv[++i];
	      Print_Grid=true;

	    }

  	  else if(temp=="-A")
	    MAX_SCORE=false;


	  else if(temp=="-m")
	    {
	      check(i,argv,argc);
	      Min_Block_Size=atoi(argv[++i]);

	    }

	  else if(temp=="-T")
	   {
	     check(i,argv,argc);
	     Threshold=atof(argv[++i]);
	     if(Threshold>1||Threshold<0)
	       {
		 cerr<<"fish:error:Threshold cannot be greater than 1 or negative\n";
                 usage();
  	       }
	   }

	   else if(temp=="-p")
	   {
	     check(i,argv,argc);
	     Threshold=atof(argv[++i]);
	     if(BLOCK_PROB>1||BLOCK_PROB<0)
	       {
		 cerr<<"fish:error:BLOCK_PROB cannot be greater than 1 or negative\n";
                 usage();
  	       }
	   }


	  else if(temp=="-S")
	   {
	     check(i,argv,argc);
	     Min_Score=atoi(argv[++i]);
	     if(Min_Score<0)
	       {
		 cerr<<"fish:error:Min_score cannot be negative\n";
		 usage();
	       }

	   }


	  else if(temp=="-H")
	   {
	     check(i,argv,argc);
	     Min_Score=atoi(argv[++i]);
	     if(Min_Score<=0)
	       {
		 cerr<<"fish:error:TOP_HITS cannot be zero or negative !\n";
		 usage();
	       }

	   }




	  else if(temp=="-D")
	  {
             check(i,argv,argc);
	     Max_Distance=atoi(argv[++i]);
	     if(Max_Distance<0)
	       {
		 cerr<<"fish:error:Max_Distance cannot be negative\n";
		 usage();
	       }


	   }

	  else
	    {
	      cerr<<"fish:error: No such option(s)\n";
	      usage();
	    }

	i++;

   }

}

void check(int j, char* argv[],int argc)
{
  if(j==argc-1)
  {
	cerr<<"fish:error:Missing arguments or wrong options!\n";
	usage();
  }

  if(*argv[++j]=='-')
  {
	cerr<<"fish:error:Missing option\n";
	usage();
   }

}

void usage(void)
{
  cerr<<"usage:  fish [-f Control_File |control.txt]\n";
  cerr<<"\t[-b Block_File] [-g Grid_File] [-C Contig_File]\n";
  cerr<<"\t[-D Max_Distance |10] [-S Min_Score |200] [-T Threshold |.05]\n";
  cerr<<"\t[-m Min_Block_Size |3] [-off Turn off Detandemize] [-h help]\n";
  cerr<<"\t[-H Top_Hits |5] [-A Average_Score Selection | Max_Score]\n";
  cerr<<"\t[-p BLOCK_PROB |.001] [-B Block_File_Simple][-q QUIET_MODE]\n";
  cerr<<"\t[-t Print timing]\n\n";
  exit(1);
}
