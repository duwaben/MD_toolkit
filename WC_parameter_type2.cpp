#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>

using namespace std;

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
	//name list
	vector<string> namelist;
	readName(namelist);
	
	NumPerType thisNumPerType;
	thisNumPerType = readNumber("XDATCAR");

	//------------------------------------------
	map<string, double> pCNinformation;
	char pCNname[64];
	double tempvalue;
	fp = fopen("pairCorrelationFunctions/pCN_information", "r");
	for(int i = 0; i < thisNumPerType.getsize(); i ++)
	{
		for(int j = i; j < thisNumPerType.getsize(); j ++)
		{
			sprintf(pCNname, "pCN_%d%d", i+1, j+1);
			findWord(pCNname);
			fscanf(fp, "%lf", &tempvalue);
			pCNinformation[pCNname] = tempvalue;
			rewind(fp);
		}
	}
	fclose(fp);

	//------------------------------------------
	map<string, double> WCinformation_type2;
	map<string, double> NormalizedWCinformation_type2;
	char WCname[64];
	double C1, C2, Z11, Z12, Z21, Z22, Z1, Z2, Zw, amax;
	for(int i = 0; i < thisNumPerType.getsize(); i ++)
	{
		for(int j = i; j < thisNumPerType.getsize(); j ++)
		{
			C1 = (double)(thisNumPerType.getmax(i)-thisNumPerType.getmin(i))/((thisNumPerType.getmax(i)-thisNumPerType.getmin(i))+(thisNumPerType.getmax(j)-thisNumPerType.getmin(j)));
			C2 = (double)(thisNumPerType.getmax(j)-thisNumPerType.getmin(j))/((thisNumPerType.getmax(i)-thisNumPerType.getmin(i))+(thisNumPerType.getmax(j)-thisNumPerType.getmin(j)));
			sprintf(pCNname, "pCN_%d%d", i+1, i+1);
			Z11 = pCNinformation[pCNname];
			sprintf(pCNname, "pCN_%d%d", i+1, j+1);
			Z12 = pCNinformation[pCNname];
			sprintf(pCNname, "pCN_%d%d", j+1, i+1);
			Z21 = pCNinformation[pCNname];
			sprintf(pCNname, "pCN_%d%d", j+1, j+1);
			Z22 = pCNinformation[pCNname];
			Z1 = Z11+Z12;
			Z2 = Z21+Z22;
			Zw = C1*Z2 + C2*Z1;
			sprintf(WCname, "a%d%d", i+1, j+1);
			WCinformation_type2[WCname] = 1-Z12/(C2*Zw);
			if(C1*Z1 > C2*Z2)
			{
				amax = 1 - Z1/(C2*Zw);
			}
			else
			{
				amax = 1 - Z2/(C1*Zw);
			}
			NormalizedWCinformation_type2[WCname] = WCinformation_type2[WCname] /amax;
		}
	}

	//output
	fp = fopen("pairCorrelationFunctions/WC_information_type2", "w");
	for(map<string, double>::iterator i = WCinformation_type2.begin(); i !=WCinformation_type2.end(); i ++)
	{
		fprintf(fp, "%s\t%lf\n", i->first.c_str(), i->second);
	}
	fclose(fp);

	fp = fopen("pairCorrelationFunctions/NormalizedWC_information_type2", "w");
	for(map<string, double>::iterator i = NormalizedWCinformation_type2.begin(); i !=NormalizedWCinformation_type2.end(); i ++)
	{
		fprintf(fp, "%s\t%lf\n", i->first.c_str(), i->second);
	}
	fclose(fp);

	return 0;
}