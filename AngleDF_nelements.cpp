#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <sys/stat.h>

using namespace std;

#define NATOM 200       //total atom number
#define NOMIT 4000 		//the number of omitted steps
#define NSTEP 8000		//the number of counting steps
#define gap 1			//the interval of counting
#define delt_a 0.5
#define MAXBIN 360
#define Pi 3.141592653

FILE *fp;
double LATTICE;
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];
double u1[NATOM], v1[NATOM], w1[NATOM];
double gr[MAXBIN];
double sumgr[MAXBIN];

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

double r(double uc1, double vc1, double wc1, double uc2, double vc2, double wc2)
{
	double xd, yd, zd;
	double rpbc;
	xd = uc1 - uc2;
	yd = vc1 - vc2;
	zd = wc1 - wc2;
	if(xd > LATTICE /2) xd -= LATTICE; if(xd < -LATTICE /2) xd += LATTICE;
	if(yd > LATTICE /2) yd -= LATTICE; if(yd < -LATTICE /2) yd += LATTICE;
	if(zd > LATTICE /2) zd -= LATTICE; if(zd < -LATTICE /2) zd += LATTICE;
	rpbc = sqrt(xd * xd + yd * yd + zd * zd);
	return rpbc;
}


double bondAngle(int atom1, int atom2, int atom3, double cutoff1, double cutoff2)
{
	double r12 = r(u1[atom1], v1[atom1], w1[atom1], u1[atom2], v1[atom2], w1[atom2]);
	double r13 = r(u1[atom1], v1[atom1], w1[atom1], u1[atom3], v1[atom3], w1[atom3]);
	if( r12 > cutoff1 || r13 > cutoff2 || r12 < 0.001 || r13 < 0.001 )
		return -1.0;
	double r23 = r(u1[atom2], v1[atom2], w1[atom2], u1[atom3], v1[atom3], w1[atom3]);
	double anglecos = ( r12 * r12 + r13 * r13 - r23 * r23 ) /( 2 * r12 * r13 );
	return acos(anglecos)*180/Pi;
}

int cmp(const pair<double, double>& x, const pair<double, double> &y)
{
	return x.second > y.second;
}

NumPerType readCoordinate(const char *name)
{
	NumPerType tempnumber;

	fp = fopen(name, "r");

	int temp, t, flag;

	const char *word1 = "Direct configuration";
	char line[512];

	//read useful informations form files
	//**********************************************************
	//**********************************************************
	fgets(line, 512, fp); fgets(line, 512, fp);
	fscanf(fp, "%lf", &LATTICE);

	for(int i = 0; i < 4; i ++)fgets(line, 512, fp);
		fgets(line, 512, fp);
	for(int i = 0; i < strlen(line); ++ i)
	{
		flag = sscanf(line+i, "%d%n", &temp, &t); i += t;
		if(flag == 1)
			tempnumber.append(temp);
	}

	//skip the steps
	for(int i = 0; i < NOMIT; i ++){
		findWord(word1);
	}
	for(int i = 0; i < NSTEP; i ++){
		findWord(word1);
		fgets(line, 512, fp);
		for(int j = 0; j < NATOM; j ++){
			//read the last 6000 steps
			fscanf(fp, "%lf%lf%lf", &u[j][i], &v[j][i], &w[j][i]);
			u[j][i] *= LATTICE;
			v[j][i] *= LATTICE;
			w[j][i] *= LATTICE;
		}
	}
	printf("Coordinate-reading part is over.\n");
	fclose(fp);
	//*********************************************************
	//*********************************************************
	return tempnumber;
}

map<string, double> readCutoff(NumPerType & temp)
{	
	char name[64];
	double cutoff;
	map<string, double> cutoffinformation;
	fp = fopen("pairCorrelationFunctions/CUTOFFPEAK_information", "r");
	for(int i = 0; i < temp.getsize(); ++ i)
	{
		for(int j = i; j < temp.getsize(); ++ j)
		{
			rewind(fp); //redirect
			sprintf(name, "CUTOFF%d%d:", i+1, j+1);
			findWord(name);
			fscanf(fp, "%lf", &cutoff);
			sprintf(name, "CUTOFF%d%d", i+1, j+1);

			cutoffinformation[name] = cutoff;
		}
	}
	fclose(fp);
	return cutoffinformation;
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


int main(int argc, char const *argv[])
{
	//namelist
	vector<string> namelist;
	readName(namelist);
	
	//
	NumPerType thisNumPerType;
	thisNumPerType = readCoordinate("XDATCAR");
	mkdir("angleDistributionFunctions", 00777);
	//cutoff
	map<string, double> cutoffinformation;
	cutoffinformation = readCutoff(thisNumPerType);


	double angle;
	int sumg_angle;
	int g_angle[MAXBIN];
	double g_angleDistribution[MAXBIN];
	char name[64];
	char filename[128];
	int BIN;

	//
	vector< pair<double, double> > angledistribution;
	int flag;
	FILE *fp1;
	fp1 = fopen("angleDistributionFunctions/SpecialList", "w");

	for(int i = 0; i < thisNumPerType.getsize(); ++ i)
	{
		for(int j = 0; j < thisNumPerType.getsize(); ++ j)
		{
			for(int k = 0; k < MAXBIN; k ++ )
			{
				g_angle[k] = 0;
			}
			sumg_angle = 0;
			angledistribution.clear();

			for(int k = 0; k < NSTEP; k += gap)
			{
				for(int l = 0; l < NATOM; ++ l)
				{
					u1[l] = u[l][k];
					v1[l] = v[l][k];
					w1[l] = w[l][k];
				}
				for(int ni = thisNumPerType.getmin(i); ni < thisNumPerType.getmax(i); ++ ni)
				{
					for(int nj = thisNumPerType.getmin(j); nj < thisNumPerType.getmax(j); ++ nj)
					{
						for(int nk = nj+1; nk < thisNumPerType.getmax(j); ++ nk)
						{
							if(i <= j)
							{
								sprintf(name, "CUTOFF%d%d", i+1, j+1);
							}
							else
							{
								sprintf(name, "CUTOFF%d%d", j+1, i+1);								
							}

							angle = bondAngle(ni, nj, nk, cutoffinformation[name], cutoffinformation[name]);
							if(angle > 0)
							{
								BIN = floor(angle /delt_a);
								g_angle[BIN] ++;
							}
						}
					}
				}

			}

			for(int k = 0; k < MAXBIN; ++ k)
			{
				sumg_angle += g_angle[k];
			}

			for(int k = 0; k < MAXBIN; ++ k)
			{
				g_angleDistribution[k] = (double)g_angle[k] /(double)sumg_angle /delt_a;
			}

			sprintf(filename, "angleDistributionFunctions/ADF%d_%d%d", i+1, j+1, j+1);
			fp = fopen(filename, "w");
			for(int k = 0; k < MAXBIN; k ++){
				fprintf(fp, "%lf\t%lf\n", (k+0.5)*delt_a, g_angleDistribution[k]);
			}			
			fclose(fp);

			//check if any highest peak appears around 95~115
			for(int k = 0; k < MAXBIN; k ++)
			{
				angledistribution.push_back(make_pair((k+0.5)*delt_a, g_angleDistribution[k]));
			}
			sort(angledistribution.begin(), angledistribution.end(), cmp);
			flag = 1;
			for(int k = 0; k < 3; k ++)
			{
				if(angledistribution[k].first < 95.0 || angledistribution[k].first > 110.0)
				{
					flag = 0;
				}
			}
			if(flag == 1)
			{
				fprintf(fp1, "%s%s%s\n", namelist[i].c_str(), namelist[j].c_str(), namelist[j].c_str());
			}
		}
	}

	fclose(fp1);
	return 0;
}
