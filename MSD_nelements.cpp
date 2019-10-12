/*
	Name: MSD(mean square displacement) Processing Program for n elements
	Date: 08/10/18 
	Input: INCAR, XDATCAR, POTCAR;
	Output: MSD_1, MSD_2, MSD_3, DMSD_information=;
	Description: This program can be used to handle ternary system
*/

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <sys/stat.h>

using namespace std;

#define TCOR 2000 //步数跨度（T0~T0+TCOR）
#define NOMIT 4000 //忽略步数
#define NSTEP 8000 //计算步数
#define NATOM 200
#define Pi  3.141593 

FILE *fp;
int TEBEG;
double LATTICE, POTIM;
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP]; 
double u1[NSTEP], v1[NSTEP], w1[NSTEP];
double MSD[TCOR], msddelta[TCOR];

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

void meanSquareDisplacement()
{ 
	int T, T0;
	double MSDtemp;
	int TMAX;
	for(T = 0; T < TCOR; T ++){
		MSDtemp = 0.0;
		TMAX = NSTEP - T;
		for(T0 = 0; T0 < TMAX; T0 ++){
			double uuu = u1[T0+T] - u1[T0];
			double vvv = v1[T0+T] - v1[T0];
			double www = w1[T0+T] - w1[T0];
						
			MSDtemp += uuu * uuu + vvv * vvv + www * www;
		}
		msddelta[T] = MSDtemp/(double)TMAX;
	}
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

void readINCAR(const char *name)
{
	const char *word1 = "POTIM";
	const char *word2 = "=";
	const char *word3 = "TEBEG";

	fp = fopen(name, "r");
	findWord(word1);
	findWord(word2);
	fscanf(fp, "%lf", &POTIM);
	fclose(fp);

	fp = fopen(name, "r");
	findWord(word3);
	findWord(word2);
	fscanf(fp, "%d", &TEBEG);
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
	readINCAR("INCAR");
	mkdir("dynamicalProperty", 00777);				


	//坐标修正：用于计算MSD和VAF
	double ru[NATOM][NSTEP], rv[NATOM][NSTEP], rw[NATOM][NSTEP];
	for(int i = 0; i < NSTEP; i ++){
		for(int j = 0; j < NATOM; j ++){
			ru[j][i] = u[j][i];
			rw[j][i] = w[j][i];
			rv[j][i] = v[j][i];
		}
	}

	for(int i = 0; i < NSTEP - 1; i ++){
		for(int j = 0; j < NATOM; j ++){
			if((ru[j][i+1] - ru[j][i]) > LATTICE/2.0)
				ru[j][i+1] -= LATTICE;
			if((ru[j][i+1] - ru[j][i]) < -LATTICE/2.0)
				ru[j][i+1] += LATTICE;
			if((rv[j][i+1] - rv[j][i]) > LATTICE/2.0)
				rv[j][i+1] -= LATTICE;
			if((rv[j][i+1] - rv[j][i]) < -LATTICE/2.0)
				rv[j][i+1] += LATTICE;
			if((rw[j][i+1] - rw[j][i]) > LATTICE/2.0)
				rw[j][i+1] -= LATTICE;
			if((rw[j][i+1] - rw[j][i]) < -LATTICE/2.0)
				rw[j][i+1] += LATTICE;
		}
	}

	char filename[128];
	char MSDname[128];
	double area, DMSD;
	double bottomEdge = TCOR * 0.001 * POTIM;
	map<string, double> MSDinformation;
	//
	for(int k = 0; k < thisNumPerType.getsize(); ++ k)
	{
		for(int T = 0; T < TCOR; T ++){
			MSD[T] = 0.0;
		}

		for(int j = thisNumPerType.getmin(k); j < thisNumPerType.getmax(k); j ++){
			for(int i = 0; i < NSTEP; i ++){
				u1[i] = ru[j][i];	
				v1[i] = rv[j][i];	
				w1[i] = rw[j][i];
			}
			meanSquareDisplacement(); 
			for(int T = 0; T < TCOR; T ++){
				MSD[T] += msddelta[T];
			}
		}

		for(int T = 0; T < TCOR; T ++){
			MSD[T] = MSD[T]/(thisNumPerType.getmax(k)-thisNumPerType.getmin(k));
		}

		sprintf(filename, "dynamicalProperty/MSD_%d", k+1);
		fp = fopen(filename, "w");
		fprintf(fp, "X\t%s\n", namelist[k].c_str());
		for(int T = 0; T < TCOR; T++){
			fprintf(fp, "%lf\t%lf\n", T * 0.001 * POTIM, MSD[T]);
		}
		fclose(fp);

		area = 0.0;
		for(int i = 0; i < TCOR; ++i)
		{
			area += 0.001 * POTIM * MSD[i];
		}
		DMSD = area *2 /(bottomEdge * bottomEdge) /6;
		sprintf(MSDname, "DMSD_%d", k+1);
		MSDinformation[MSDname] = DMSD;
	}

	fp = fopen("dynamicalProperty/DMSD_information", "w");
	for(map<string, double>::iterator i = MSDinformation.begin(); i != MSDinformation.end(); ++ i)
	{
		fprintf(fp, "%s\t%lf\n", (*i).first.c_str(), (*i).second);
	}
	fclose(fp);
	
	return 0;	
}
