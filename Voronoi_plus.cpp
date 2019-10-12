/*
	The abundance and lifetime of Voronoi polyhedra


*/

#include <iostream>
#include <vector>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>


using namespace std;

#define NSTEP_voro 3000
#define NATOM 200

FILE *fp;

class NumPerType
{
public:
	NumPerType();
	int getmin(int);	
	int getmax(int);
	int total();
	void append(int);
	int getsize();
private:
	vector<int> number;
};

NumPerType::NumPerType()
{
}

int NumPerType::getmin(int type)
{
	int sum = 0;
	for(int i = 0; i < type; ++ i)
	{
		sum += number[i];
	}
	return sum;
}

int NumPerType::getmax(int type)
{
	int sum = 0;
	for(int i = 0; i < type+1; ++ i)
	{
		sum += number[i];
	}
	return sum;
	
}

int NumPerType::total()
{
	int sum = 0;
	for(int i = 0; i < number.size(); ++ i)
	{
		sum += number[i];
	}
	return sum;

}

void NumPerType::append(int temp)
{
	number.push_back(temp);
}

int NumPerType::getsize()
{
	return number.size();
}

void findWord(const char *word)
{
	int n = strlen(word);
	char word0[n];
	int i;
	for(i = 0; i < n; i ++)
		fscanf(fp, "%c", &word0[i]);
	while(strncmp(word0, word, n) != 0){
		for(i = 0; i < n-1; i ++)
			word0[i] = word0[i+1];
		fscanf(fp, "%c", &word0[n-1]);
	}	
}


NumPerType readNumber(const char *name)
{
	NumPerType tempnumber;

	fp = fopen(name, "r");

	int temp, t, flag;

	const char *word1 = "Direct configuration";
	char line[512];

	//read useful informations form files
	//**********************************************************
	//**********************************************************
	
	for(int i = 0; i < 7; i ++)fgets(line, 512, fp);
	for(int i = 0; i < strlen(line); ++ i)
	{
		flag = sscanf(line+i, "%d%n", &temp, &t); i += t;
		if(flag == 1)
			tempnumber.append(temp);
	}

	fclose(fp);
	//*********************************************************
	//*********************************************************
	return tempnumber;
}

void readName(vector<string> & namelist)
{
	fp = fopen("XDATCAR", "r");
	char line[256];
	char name[256];
	int t, flag;
	string namestring;
	for(int i = 0; i < 6; ++ i)
	{
		fgets(line, 256, fp);
	}
	for(int i = 0; i < strlen(line); ++ i)
	{
		flag = sscanf(line+i, "%s%n", name, &t); i += t;
		namestring = name;
		if(flag == 1)
			namelist.push_back(namestring);
	}
	fclose(fp);
}

void getLifeTime(map< string, pair<int, int> > &lifeTimeCount)
{
	const char *word = "Freq_table =";
	string index_pre = "0";
	string index;
	int freq_table[7];
	char buffer[128];
	fp = fopen("voronoiTessellation/Voronoi_EachAtom", "r");
	for(int i = 0; i < NSTEP_voro; i ++)
	{
		findWord(word);
		for(int j = 0; j < 6; j ++)
		{
			fscanf(fp, "%d", &freq_table[j]);
		}
		if(fscanf(fp, "%d", &freq_table[6]) == 0)freq_table[6] = 0;

		sprintf(buffer, "<%d,%d,%d,%d>", freq_table[3],freq_table[4],freq_table[5],freq_table[6]);
		index = buffer;

		if(index != index_pre)
		{
			if(lifeTimeCount.find(index) == lifeTimeCount.end())
			{
				lifeTimeCount[index].first = 1;
				lifeTimeCount[index].second = 1;				
			}
			else
			{
				lifeTimeCount[index].first ++;
				lifeTimeCount[index].second ++;
			}
		}
		else
		{
			lifeTimeCount[index].first ++;
		}

		index_pre = index;	
	}
	fclose(fp);
}

int cmp(const pair<string, double>& x, const pair<string, double>& y)    
{    
    return x.second > y.second;    
}

int main(int argc, char const *argv[])
{
	//namelist
	vector<string> namelist;
	readName(namelist);

	NumPerType thisNumPerType;
	thisNumPerType = readNumber("XDATCAR");

	char line[1024];
	char filename[128];
	int Atom_id;
	map< string, pair<int, int> > lifeTimeCount;
	map< string, pair<int, int> > lifeTimeCount_total;
	vector< pair<string, double> > lifetime;
	vector< pair<string, double> > lifetime_total;  
	FILE *fp1; FILE *fp2;
	vector<double> averageLifeTime;
	int sumtime, Voronoinumber;
	//*************************************************************
	//*************************************************************
	fp1 = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	for(int i = 0; i < thisNumPerType.getsize(); i ++)
	{
		for(int j = thisNumPerType.getmin(i); j < thisNumPerType.getmax(i); j ++)
		{
			fp2 = fopen("voronoiTessellation/Voronoi_EachAtom", "w");
			while(fgets(line, 1024, fp1))
			{
				sscanf(line, "%*[^=]=%d", &Atom_id);
				if(Atom_id == j+1)
				{
					fputs(line, fp2);
				}
			}
			fclose(fp2);
			rewind(fp1);

			getLifeTime(lifeTimeCount);
			getLifeTime(lifeTimeCount_total);
		}

		//get average value
		sumtime = Voronoinumber = 0;
		for(map< string, pair<int, int> >::iterator iter = lifeTimeCount.begin(); iter != lifeTimeCount.end(); iter ++)
		{
			sumtime += iter->second.first;
			Voronoinumber += iter->second.second;			
		}
		averageLifeTime.push_back((double)sumtime/Voronoinumber);

		//sort
		for(map< string, pair<int, int> >::iterator iter = lifeTimeCount.begin(); iter != lifeTimeCount.end(); iter ++)
		{
			lifetime.push_back(make_pair(iter->first, (double)iter->second.first/iter->second.second));
		}
		sort(lifetime.begin(), lifetime.end(), cmp);

		//output
		//-------------------------------------------------------------------------------------------------------------
		sprintf(filename, "voronoiTessellation/Lifetime_%d_%s", i+1, namelist[i].c_str());
		fp2 = fopen(filename, "w");
		fprintf(fp2, "X\t%s\n", namelist[i].c_str());
		for(int j = 0; j < lifetime.size(); ++j)
		{
			fprintf(fp2, "%s\t%lf\n", lifetime[j].first.c_str(), lifetime[j].second);
		}
		fclose(fp2);
		//-------------------------------------------------------------------------------------------------------------
		lifeTimeCount.clear();
		lifetime.clear();
	}

	fclose(fp1);

	//get total average value
	sumtime = Voronoinumber = 0;
	for(map< string, pair<int, int> >::iterator iter = lifeTimeCount_total.begin(); iter != lifeTimeCount_total.end(); iter ++)
	{
		sumtime += iter->second.first;
		Voronoinumber += iter->second.second;			
	}
	averageLifeTime.push_back((double)sumtime/Voronoinumber);

	//sort
	for(map< string, pair<int, int> >::iterator iter = lifeTimeCount_total.begin(); iter != lifeTimeCount_total.end(); iter ++)
	{
		lifetime_total.push_back(make_pair(iter->first, (double)iter->second.first/iter->second.second));
	}
	sort(lifetime_total.begin(), lifetime_total.end(), cmp);

	//output
	//-------------------------------------------------------------------------------------------------------------
	fp2 = fopen("voronoiTessellation/Lifetime_total", "w");
	fprintf(fp2, "X\ttotal\n");
	for(int j = 0; j < lifetime_total.size(); ++j)
	{
		fprintf(fp2, "%s\t%lf\n", lifetime_total[j].first.c_str(), lifetime_total[j].second);
	}
	fclose(fp2);

	//output average
	fp2 = fopen("voronoiTessellation/Lifetime_average", "w");
	for(int i = 0; i < averageLifeTime.size()-1; i ++)
	{
		fprintf(fp2, "%s\t%lf\n", namelist[i].c_str(),averageLifeTime[i]);
	}
	fprintf(fp2, "total\t%lf\n", averageLifeTime[averageLifeTime.size()-1]);
	fclose(fp2);

	return 0;
}
