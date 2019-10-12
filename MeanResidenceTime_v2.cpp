#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#define NSTEP_voro 3000
#define NATOM 200

FILE *fp;
double POTIM;

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

//functional function
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

void readINCAR()
{
	char *word1 = "POTIM";
	char *word2 = "=";

	fp = fopen("INCAR", "r");
	findWord(word1);
	findWord(word2);
	fscanf(fp, "%lf", &POTIM);
	fclose(fp);
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

void getLifeTime_bs(map< string, pair<int, int> > &lifeTimeCount, NumPerType thisNumPerType, int atomType1, int atomType2, char* bondName)
{
	vector<int> index_pre;
	vector<int> index;
	char buffer[256];
	char neighbourstring[256];
	char line[1024];
	vector<int> tempnumber;
	int flag, temp, t;
	int newBondNum;
	fp = fopen("voronoiTessellation/Voronoi_EachAtom", "r");
	for(int i = 0; i < NSTEP_voro; i ++)
	{
		fgets(line, 1024, fp);
		sscanf(line, "%*[^=]=%*[^=]=%*[^=]=%[^|]", neighbourstring);
		tempnumber.clear();

		for(int k = 0; k < strlen(neighbourstring); k ++)
		{
			flag = sscanf(neighbourstring+k, "%d%n", &temp, &t); k+=t;
			if(flag == 1 && temp-1 >= thisNumPerType.getmin(atomType2) && temp-1 < thisNumPerType.getmax(atomType2))
				tempnumber.push_back(temp);
		}

		sort(tempnumber.begin(), tempnumber.end());

		index = tempnumber;
		//minus
		newBondNum = 0;
		for(int k = 0; k < index.size(); k ++)
		{
			flag = 0;
			for(int l = 0; l < index_pre.size() && flag == 0; l ++)
			{
				if(index[k] == index_pre[l])
				{
					flag = 1;
				}
			}
			if(flag == 0){
				newBondNum ++;
			}
		}

		if(lifeTimeCount.find(bondName) == lifeTimeCount.end())
		{
			lifeTimeCount[bondName].first = index.size();
			lifeTimeCount[bondName].second = newBondNum;
		}
		else
		{
			lifeTimeCount[bondName].first += index.size();
			lifeTimeCount[bondName].second += newBondNum;
		}
		
		index_pre = index;
	}
	fclose(fp);
}

int main(int argc, char const *argv[])
{
	//name list
	vector<string> namelist;
	readName(namelist);
	
	NumPerType thisNumPerType;
	thisNumPerType = readNumber("XDATCAR");

	readINCAR();
	//---------------------------------------------------------

	//lifetime of bonds
	FILE *fp1; FILE *fp2;
	map< string, pair<int, int> > lifeTimeCount;
	vector< pair<string, double> > lifetime;
	int Atom_id;
	char bondName[128];
	char line[1024];
	//--------------------------------------------------------------------------------------
	fp1 = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	for(int i = 0; i < thisNumPerType.getsize(); i ++)
	{
		for(int j = i; j < thisNumPerType.getsize(); j ++)
		{
			sprintf(bondName, "%s-%s", namelist[i].c_str(), namelist[j].c_str());

			for(int k = thisNumPerType.getmin(i); k < thisNumPerType.getmax(i); k ++)
			{
				fp2 = fopen("voronoiTessellation/Voronoi_EachAtom", "w");
				while(fgets(line, 1024, fp1))
				{
					sscanf(line, "%*[^=]=%d", &Atom_id);
					if(Atom_id == k+1)
					{
						fputs(line, fp2);
					}
				}
				fclose(fp2);
				rewind(fp1);

				getLifeTime_bs(lifeTimeCount, thisNumPerType, i, j, bondName);
			}
		}
	}

	//
	for(map< string, pair<int, int> >::iterator iter = lifeTimeCount.begin(); iter != lifeTimeCount.end(); iter ++)
	{
		lifetime.push_back(make_pair(iter->first, (double)iter->second.first/iter->second.second));
	}
	fclose(fp1);

	//output
	fp = fopen("voronoiTessellation/BondStrength", "w");
	for(int i = 0; i < thisNumPerType.getsize(); i ++)
	{
		for(int j = i; j < thisNumPerType.getsize(); j ++)
		{
			sprintf(bondName, "%s-%s", namelist[i].c_str(), namelist[j].c_str());
			for(int k = 0; k < lifetime.size(); k ++)
			{
				if(lifetime[k].first == bondName)
				{
					fprintf(fp, "%s\t%lf\n", bondName, lifetime[k].second * POTIM);
				}
			}
		}
	}
	fclose(fp);

	return 0;
}