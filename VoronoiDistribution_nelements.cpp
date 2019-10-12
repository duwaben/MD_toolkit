#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>

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

int cmp(const pair<string, int>& x, const pair<string, int>& y)    
{    
    return x.second > y.second;    
}    
     
void sortMapByValue(map<string, int>& tMap, vector< pair<string, int> >& tVector)    
{    
    for (map<string, int>::iterator curr = tMap.begin(); curr != tMap.end(); curr++)     
        tVector.push_back(make_pair(curr->first, curr->second));      
     
    sort(tVector.begin(), tVector.end(), cmp);    
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

int main(int argc, char const *argv[])
{
	//namelist
	vector<string> namelist;
	readName(namelist);

	NumPerType thisNumPerType;
	thisNumPerType = readNumber("XDATCAR");
	
	const char *word1 = "ID =";
	const char *word2 = "Freq_table =";
	int freq_table[7];
	char buffer[128];

	map<string, int> dic;
	vector< pair<string, int> > tVector;    
	//****************************************************************************************
	//****************************************************************************************
	fp = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	for(int i = 0; i < NSTEP_voro*NATOM; ++ i){
		findWord(word2);
		for(int j = 0; j < 6; ++ j){
			fscanf(fp, "%d", &freq_table[j]);
		}

		if(fscanf(fp, "%d", &freq_table[6]) == 0)freq_table[6] = 0;

		sprintf(buffer, "<%d,%d,%d,%d>", freq_table[3],freq_table[4],freq_table[5],freq_table[6]);
		if(dic.find(buffer) == dic.end())
			dic[buffer] = 1;
		else 
			dic[buffer] ++;
	}
	fclose(fp);

	//sort
    sortMapByValue(dic, tVector);
	
	fp = fopen("voronoiTessellation/VoronoiDistribution_total_information", "w");
	for(int i = 0; i < tVector.size(); i ++)  
		fprintf(fp, "%s || %lf || %d\n", tVector[i].first.c_str(), double(tVector[i].second)/(NSTEP_voro*NATOM), tVector[i].second);
	fclose(fp);

	//-----------------------------------------------------------------------------------
	int v1, v2, v3, v4;
	int number, sumnumber;
	int hist[21];
	char line[128];

	for(int i = 0; i < 21; ++ i)
	{
		hist[i] = 0;
	}
	sumnumber = 0;

	fp = fopen("voronoiTessellation/VoronoiDistribution_total_information", "r");
	while(fgets(line, 128, fp)){
		sscanf(line, "<%d,%d,%d,%d>||%*[^|]||%d", &v1, &v2, &v3, &v4, &number);		
		hist[v1+v2+v3+v4] += number;
		sumnumber += number;		
	}
	fclose(fp);

	fp = fopen("voronoiTessellation/CNDistribution_total_information", "w");
	for(int i = 0; i < 21; ++ i)
		fprintf(fp, "%d\t%lf\n", i, (double)hist[i]/sumnumber);	
	fclose(fp);
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	int atom_id, totalnumber;
	char filename[128];
	for(int k = 0; k < thisNumPerType.getsize(); ++ k)
	{
		dic.clear();
		tVector.clear();
		fp = fopen("voronoiTessellation/VORONOI_RESULT", "r");
		totalnumber = 0;
		for(int i = 0; i < NSTEP_voro*NATOM; ++ i){
			
			findWord(word1);
			fscanf(fp, "%d", &atom_id);
			if(atom_id > thisNumPerType.getmin(k) && atom_id <= thisNumPerType.getmax(k)){
				totalnumber ++;
				findWord(word2);
				for(int j = 0; j < 6; ++ j){
					fscanf(fp, "%d", &freq_table[j]);
				}

				if(fscanf(fp, "%d", &freq_table[6]) == 0)freq_table[6] = 0;

				sprintf(buffer, "<%d,%d,%d,%d>", freq_table[3],freq_table[4],freq_table[5],freq_table[6]);
				if(dic.find(buffer) == dic.end())
					dic[buffer] = 1;
				else 
					dic[buffer] ++;
			}
		}
		fclose(fp);
 		
 		//sort  
    	sortMapByValue(dic, tVector);
    	sprintf(filename, "voronoiTessellation/VoronoiDistribution_%d_%s_information", k+1, namelist[k].c_str());
    	fp = fopen(filename, "w");
		for(int i = 0; i < tVector.size(); i ++)  
			fprintf(fp, "%s || %lf || %d\n", tVector[i].first.c_str(), double(tVector[i].second)/totalnumber, tVector[i].second);
		fclose(fp);

		//CNdistribution
		for(int i = 0; i < 21; ++ i)
		{
			hist[i] = 0;
		}
		sumnumber = 0;

		fp = fopen(filename, "r");
		while(fgets(line, 128, fp)){
			sscanf(line, "<%d,%d,%d,%d>||%*[^|]||%d", &v1, &v2, &v3, &v4, &number);		
			hist[v1+v2+v3+v4] += number;
			sumnumber += number;		
		}
		fclose(fp);

		sprintf(filename, "voronoiTessellation/CNDistribution_%d_%s_information", k+1, namelist[k].c_str());
		fp = fopen(filename, "w");
    	fprintf(fp, "X\t%s\n", namelist[k].c_str());
		for(int i = 0; i < 21; ++ i)
			fprintf(fp, "%d\t%lf\n", i, (double)hist[i]/sumnumber);	
		fclose(fp);		
	}

	return 0;
}
