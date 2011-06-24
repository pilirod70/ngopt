void dpalgo(void)
{

int i,j;
double l;

typedef multiset<TL1,sort_crit> SET1; //multiset declared for storing highest scoring nodes

SET1 E_P; //set of endpoints with higest scores


/*forward steps of the dp algorithm*/

for(i=0;i<EX.size();i++)
{
	if(EX[i].DP==0)
		EX[i].DP=1;

	for(j=i+1;j<EX.size();j++)
	{
		if(EX[j].feat.y-EX[i].feat.y>Thresh_Dist-1)
		break;



		if((EX[i].feat.y<EX[j].feat.y)&&(EX[i].feat.x!=EX[j].feat.x)&&((abs(EX[i].feat.x-EX[j].feat.x)+EX[j].feat.y-EX[i].feat.y)<=Thresh_Dist))

		{

			if(EX[j].DP==0)
			{
				EX[j].DP=EX[i].DP +1;
				EX[j].NH.insert(EX[j].NH.end(),i);
				EX[i].in=EX[i].in+1;
				//EX[j].out=EX[j].NH.size();
				continue;
			}
			else if(EX[j].DP<EX[i].DP+1)
			{

				EX[j].NH.clear();
				EX[j].NH.insert(EX[j].NH.end(),i);
				EX[j].DP=EX[i].DP+1;
				EX[i].in=EX[i].in+1;
				continue;

			}

			else if(EX[j].DP==EX[i].DP+1)//the preference rules apply here
			{
				//EX[j].NH.insert(EX[j].NH.end(),i);
				//EX[j].out=EX[j].NH.size();
				continue;
			}

		}

	}

}

BLOCKS temp;
TL1 temp1;
list<int>::iterator pos;
list<int>::iterator pos1;
int idx;


/*loop for collecting all end points with decreasing order of score*/

	for(j=EX.size()-1;j>=0;j--)
	{
		if((EX[j].NH.size()>0)&&(EX[j].in==0))
		{
			temp1.score=EX[j].DP;
			temp1.idx=j;
			E_P.insert(temp1);
		}
	}


SET1::iterator pos3;

/*trace back of blocks of highest score*/

	for(pos3=E_P.begin();pos3!=E_P.end();++pos3)
	{
		idx=pos3->idx;

		temp.block.insert(temp.block.end(),EX[idx].index);

		pos=EX[idx].NH.begin();
		idx=*pos;


		while(1)
		{
			if(EX[idx].feat.check==false)
			{
				temp.block.insert(temp.block.end(),EX[idx].index);
				EX[idx].feat.check=true;

				if(EX[idx].NH.size()==0)
					break;

				pos=EX[idx].NH.begin();
				idx=*pos;
			}
			else
				break;

		}

		if(temp.block.size()>=Min_Block_Size)
		{
			MAX_BLOCK=(MAX(MAX_BLOCK,temp.block.size()));
			temp.flag=false;
			temp.block.reverse();
			block_list.push_back(temp);
			block_counter++;
		}


		temp.block.clear();

	}



temp.block.clear();


E_P.clear();

}//end of function dpalgo.



