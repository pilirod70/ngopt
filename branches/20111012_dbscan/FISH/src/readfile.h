/*
the multimap will contain hash keys which is
a string value like 20x3 etc. The sort criterion of the multimap
is based on a string comparison.
Specifications for the control file:
1. the contig1 vs contig2 comparison must have a file
in which contig1<= contig2 or else no comparison is made even if
the file contig2 vs contig1 is present */



/*Declaration of Functions for reading control file*/

void break_string(string s1,string& s2,string& s3); /*Function for parsing the string*/
string i2s(int x);                                  /*Function for int to string*/
void read_controlfile(void)
{
 char Tgene[MAXBUFFER];
 string gene_name;
 char xval[MAXBUFFER], yval[MAXBUFFER];   /*variable for x and y coordinates for the gene names*/
 ID2N temp_id;
 float pvalue;
 CHROMO temp;
 point temp_point;
 int max_feature;
 FILE *file_name;
 string s;            /*string variable to store lines from files*/
 ifstream C_FILE;

    C_FILE.open(Control_File.c_str(),ios::in);


 if(C_FILE.fail())
        {
   cerr<<"No Control File Specified\n";
          exit(1);
        }
        
     
 genes Ngene;       /*Ngene is of type genes for storing direction and index*/
 string contig_number;
 GENE_LIST gene_list;
 while(1)
 {
  getline(C_FILE,s);
  
  if(s=="-maps")
   continue;//getline(C_FILE,s);

  if(s=="-matches")
   break;
  else
   {
    int gene,direction,c;
    istringstream instr(s);
    string path,str;
    //string j="joy";
    instr>>contig_number>>path;
    file_name=fopen(path.c_str(),"r");
    if(file_name == NULL)
           {
       cerr<<"Map File '"<<path.c_str()<<"' not found. Please check the control file.\n";
              exit(1);
           }
    //printf("\n contig %s  path  \n",path.c_str());

    max_feature= 0;
    while(1)
    {
    fscanf(file_name,"%s %d",&Tgene,&Ngene.direction);

     if((c=getc(file_name))==EOF)
      break;

     max_feature++;
     gene_name=string(Tgene);
     temp_id.gene_name=string(Tgene);
     temp_id.direction=Ngene.direction;
     Ngene.idx=max_feature;
     gene_name=contig_number+gene_name;
     gene_list.insert(make_pair(gene_name,Ngene));
     id2n.push_back(temp_id);
    }
    fclose(file_name);


   }
   temp.index=atoi(contig_number.c_str());
   temp.feature=0;
   temp.MAX_N=max_feature;
   MaxGene=MaxGene+max_feature;
   nMAX.push_back(temp);

 }


GENE_LIST::iterator pos12;

GENE_LIST::iterator chk_pos1;
GENE_LIST::iterator chk_pos2;



 int counter=0;
 nSHIFT.push_back(0);
 for(int i=0;i<nMAX.size()-1;i++)
 {
  counter=nMAX[i].MAX_N+counter;
  nSHIFT.push_back(counter);

 }


 string Path_Name;
 int inserted=0;   //to check whether any matches have been actually inserted in the GENE
                   //data structure
 FILE *fp;
 SUBGRID temp_subgrid;


 /*Begin reading of match files*/


 SYMM temp_map;    /*temporary map for storing image match files*/
 SYMM::iterator pos2;
 MAP_LIST image_maps;
 MAP_LIST::iterator pos1;


 while(1)
 {
  getline(C_FILE,s);

  if(s=="end")
   break;

  istringstream CONTIG(s);
  string contig1,contig2;
  CONTIG>>contig1>>contig2>>Path_Name;
  int score=0;
  inserted=0;
  int int_contig1,int_contig2;
  string key_x,key_y;
  int_contig1=atoi(contig1.c_str());
  int_contig2=atoi(contig2.c_str());
  int tempscore=0;

  /*This part of the programs reads self matches of type 1vs.1 or 2vs.2 etc.
  and takes care of symmetry as well.*/

  if(int_contig1==int_contig2)
  {
   temp_subgrid.contig1=int_contig1;
   temp_subgrid.contig2=int_contig2;
   temp_subgrid.ext=false;
   temp_subgrid.points=0;
   sub_grid.push_back(temp_subgrid);

   fp=fopen(Path_Name.c_str(),"r");
   if(fp == NULL)
          {
      cerr<<"Match File '"<<Path_Name.c_str()<<"' not found.  Please check the control file.\n";
             exit(1);
          }

   while(1)
   { int c;

    fscanf(fp,"%s %s %d",&yval,&xval,&score);


    if((c=getc(fp))==EOF)
    {
      if(temp_map.size()==0)
      cout<<"Warning! No matches were included from comparison of "<<contig1<<" vs. "<<contig2<<endl;

      break;
    }

   key_y=string(yval);
   key_x=string(xval);

   //Make sure that a map entry exists before inserting the match into anything
   chk_pos1 = gene_list.find(contig1+key_y);
   chk_pos2 = gene_list.find(contig1+key_x);
   if(chk_pos1 == gene_list.end()) 
   {
      cout<<"Warning! No no map entry for:"<<key_y<<" contig id:  "<<contig1<<endl;
   }
   if(chk_pos2 == gene_list.end())
   {
      cout<<"Warning! No no map entry for:"<<key_x<<" contig id:  "<<contig1<<endl;
   }

    if((score>=Min_Score) && 
       (chk_pos1 != gene_list.end()) &&  
       (chk_pos2 != gene_list.end()))
    {

     //if(atoi(yval)<atoi(xval))
     if(gene_list[contig1+key_y].idx<gene_list[contig1+key_x].idx)
     {
        inserted++;
        tempscore=temp_map[contig1+key_y+'&'+contig2+key_x];

        if(MAX_SCORE)
          temp_map[contig1+key_y+'&'+contig2+key_x]=(MAX(score,tempscore));
        else
        {
          if(tempscore==0)
            temp_map[contig1+key_y+'&'+contig2+key_x]=(score);
          else
            temp_map[contig1+key_y+'&'+contig2+key_x]=(AVG(tempscore,score));
        }
     }

     //if(atoi(yval)>atoi(xval))
     if(gene_list[contig1+key_y].idx>gene_list[contig1+key_x].idx)
     {
        inserted++;
        tempscore=temp_map[contig1+key_x+'&'+contig2+key_y];

        if(MAX_SCORE)
          temp_map[contig1+key_x+'&'+contig2+key_y]=(MAX(score,tempscore));
        else
        {
             if(tempscore==0)
               temp_map[contig1+key_x+'&'+contig2+key_y]=(score);
             else
               temp_map[contig1+key_x+'&'+contig2+key_y]=(AVG(tempscore,score));
        }
     }

    }



   }

   if(temp_map.size()>0)
   {
     sub_grid[sub_grid.size()-1].ext=true;
     image_maps.insert(make_pair(contig1+"&"+contig2,temp_map));
   }

   temp_map.clear();
   fclose(fp);


  }



  bool mapexists=false;
  bool FIRST=false, SECOND=false;


  /*this part reads matches of the form 1vs.2 2vs.1 etc.*/
  /*IT IS NECESSARY TO HAVE AT LEAST ONE COMPARISON OF THE FORM Chromosom A VS. Chromosome B IN THE CONTROL FILE
  WHERE A<=B OTHERWISE IT WON'T READ AT ALL*/

  if((FIRST=int_contig1<int_contig2)||(SECOND=int_contig1>int_contig2))
  {
   if(FIRST)
   {
   temp_subgrid.contig1=int_contig1;
   temp_subgrid.contig2=int_contig2;
   temp_subgrid.ext=false;
   temp_subgrid.points=0;
   sub_grid.push_back(temp_subgrid);
   }
   if(SECOND)
   {
     swap(contig1,contig2);
   }



   fp=fopen(Path_Name.c_str(),"r");
   if(fp == NULL)
          {
      cerr<<"Match File '"<<Path_Name.c_str()<<"' not found.  Please check the control file.\n";
             exit(1);
          }
    while(1)
    {
     int c;
     fscanf(fp,"%s %s %d",&yval,&xval,&score);

     if((c=getc(fp))==EOF)
     {
      if(inserted==0)
      cout<<"Warning! No matches were included from comparison of: "<<contig1<<" vs. "<<contig2<<endl;

      break;
     }


     //verifiy a map entry for both matches
     key_y=string(yval);
     key_x=string(xval);
     if(SECOND)
     {
       swap(key_y,key_x);
     }


     chk_pos1 = gene_list.find(contig1+key_y);
     chk_pos2 = gene_list.find(contig2+key_x);
     if(chk_pos1 == gene_list.end()) 
     {
        cout<<"Warning! No no map entry for:"<<key_y<<" contig id:  "<<contig1<<endl;
     }
     if(chk_pos2 == gene_list.end())
     {
        cout<<"Warning! No no map entry for:"<<key_x<<" contig id:  "<<contig2<<endl;
     }

    if((score>=Min_Score) && 
       ( chk_pos1 != gene_list.end()) &&  
       ( chk_pos2 != gene_list.end()))
     {

   //     if(SECOND)
    //      swap(key_y,key_x);

        if(image_maps.size()>0)
        {

          pos1=image_maps.find(contig1+"&"+contig2);

          if(pos1==image_maps.end())
          {
           inserted++;
           mapexists=false;
           temp_map.insert(make_pair(contig1+key_y+'&'+contig2+key_x,score));
          }
          else
          {
            inserted++;
            mapexists=true;
            tempscore=(image_maps[contig1+"&"+contig2][contig1+key_y+'&'+contig2+key_x]);

            if(MAX_SCORE==true)
              image_maps[contig1+"&"+contig2][contig1+key_y+'&'+contig2+key_x]=(MAX(tempscore,score));
            else
            {
               if(tempscore==0)
                 image_maps[contig1+"&"+contig2][contig1+key_y+'&'+contig2+key_x]=(score);
               else
                 image_maps[contig1+"&"+contig2][contig1+key_y+'&'+contig2+key_x]=(AVG(tempscore,score));
            }

           }

          }

          else
          {
            temp_map.insert(make_pair(contig1+key_y+'&'+contig2+key_x,score));
            inserted++;
          }
     }

    }

    if(mapexists==false)
    image_maps.insert(make_pair(contig1+"&"+contig2,temp_map));

    temp_map.clear();

    fclose(fp);

    if(inserted>0&&sub_grid.size()>0&&FIRST)
    {
      sub_grid[sub_grid.size()-1].ext=true;

    }

  }

 }

 C_FILE.close();



 sort(sub_grid.begin(),sub_grid.end(),subgrid_sort);


 /*Copying the matches to the GENE data structure*/

 vector<point>::iterator pos3;

 string contig1,contig2,contigs;
 string key1,key2;
 for(int k=0;k<sub_grid.size();k++)
 {
  if(sub_grid[k].ext)
  {
   contig1=i2s(sub_grid[k].contig1);
   contig2=i2s(sub_grid[k].contig2);
   contigs=contig1+"&"+contig2;
   for(pos2=image_maps[contigs].begin();pos2!=image_maps[contigs].end();++pos2)
   {
   //cout<<"The firs key "<<pos2->first;
   break_string(pos2->first, key1, key2);
   //cout<<"  the corresponding keys are "<<key1<<"  and   "<<key2<<endl;
   temp_point.y=gene_list[key1].idx+nSHIFT[sub_grid[k].contig1-1];
   temp_point.x=gene_list[key2].idx+nSHIFT[sub_grid[k].contig2-1];
   temp_point.mark=false;
   temp_point.match_score=pos2->second;
   temp_point.direction=(gene_list[key1].direction)*(gene_list[key2].direction);
   GENE.push_back(temp_point);
   INV_GENE.push_back(temp_point);
   }

   sort(INV_GENE.begin(),INV_GENE.end(),sortcriterion);
   CG.insert(CG.end(),INV_GENE.begin(),INV_GENE.end());
   INV_GENE.clear();
   CG[CG.size()-1].mark=true;
   //GENE[GENE.size()-1].mark=true;


  }
 }

if(GENE.size()==0 ||GENE.size()<=2)
 {
   cerr<<"\nWarning: no matches or too few matches were included from match files please refer manual.\n";
   exit(1);
 }



//image_maps.~map();

 INV_GENE.clear();
 INV_GENE=GENE;
 //CG=GENE;   /*Copying the contents of GENE to CG for preserving the order
   //         in which matches are read*/






}


void break_string(string s, string& key1, string& key2)
{ int pos;
 int length=s.length();
 for(pos=0;pos<length;pos++)
 if(s[pos]=='&')
 break;

 key1=s.substr(0,pos);
 key2=s.substr(pos+1);
}

string i2s(int x)
{
   stringstream ss;
   ss<<x;
   string s;
   ss>>s;
   return s;
}



