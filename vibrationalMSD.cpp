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

#define TCOR 21 //步数跨度（T0~T0+20）
#define NOMIT 0 //忽略步数
#define NSTEP 80 //计算步数
#define NATOM 200
#define Pi  3.141593 

FILE *fp;
int TEBEG;
double LATTICE, POTIM;
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP]; 
double ru[NATOM][NSTEP], rv[NATOM][NSTEP], rw[NATOM][NSTEP];

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

NumPerType getNumber(const char *name)
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
	fclose(fp);
	//*********************************************************
	//*********************************************************
	return tempnumber;
}

void readName(vector<string> & namelist)
{
	fp = fopen("POSCAR", "r");
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

void readCoordination(char *filename)
{
	fp = fopen(filename, "r");
	char line[256];
	const char *word = "Direct configuration";
	
	for(int i = 0; i < NOMIT; ++ i)
	{
		findWord(word);
	}

	for(int i = 0; i < NSTEP; ++ i)
	{
		findWord(word);
		fgets(line, 256, fp);
		for(int j = 0; j < NATOM; ++ j)
		{
			fscanf(fp, "%lf%lf%lf", &u[j][i], &v[j][i], &w[j][i]);
			u[j][i] *= LATTICE;
			v[j][i] *= LATTICE;
			w[j][i] *= LATTICE;
		}
	}

	fclose(fp);
}

void modifyCoordination()
{
	for(int i = 0; i < NSTEP; i ++){
		for(int j = 0; j < NATOM; j ++){
			ru[j][i] = u[j][i];
			rv[j][i] = v[j][i];
			rw[j][i] = w[j][i];
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
}

double vibrationalmsd(vector<double> &perrun, int atomType, NumPerType thisNumPerType)
{
	vector<double> msdPerAtom;
	double temp;
	double totalPerAtom, totalPerRun, averagePerRun;
	double ruEqui[NATOM], rvEqui[NATOM], rwEqui[NATOM];
	for(int i = 0; i < NATOM; ++ i)
	{
		ruEqui[i] = rvEqui[i] = rwEqui[i] = 0;		
	}

	for(int i = 0; i < NATOM; ++ i)
	{
		for(int j = 0; j < TCOR; ++ j)
		{
			ruEqui[i] += ru[i][j];
			rvEqui[i] += rv[i][j];
			rwEqui[i] += rw[i][j];
		}
		ruEqui[i] = ruEqui[i]/TCOR; 
		rvEqui[i] = rvEqui[i]/TCOR; 
		rwEqui[i] = rwEqui[i]/TCOR; 
	}

	for(int i = thisNumPerType.getmin(atomType); i < thisNumPerType.getmax(atomType); ++ i)
	{
		for(int j = 0; j < TCOR; ++ j)
		{
			temp = (ru[i][j]-ruEqui[i])*(ru[i][j]-ruEqui[i]) + (rv[i][j]-rvEqui[i])*(rv[i][j]-rvEqui[i]) + (rw[i][j]-rwEqui[i])*(rw[i][j]-rwEqui[i]);
			msdPerAtom.push_back(temp);
		}

		totalPerAtom = 0;
		for(int j = 0; j < msdPerAtom.size(); ++ j)
		{
			totalPerAtom += msdPerAtom[j];
		}
		temp = totalPerAtom /msdPerAtom.size();
		perrun.push_back(temp);

		msdPerAtom.clear();
	}

	for(int i = 0; i < perrun.size(); ++ i)
	{
		totalPerRun += perrun[i];
	}
	averagePerRun = totalPerRun /perrun.size();

	//average over different atoms
	return averagePerRun;
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

	//------------------------------------------
	NumPerType thisNumPerType;
	thisNumPerType = getNumber("POSCAR");
	readINCAR("INCAR");

	//------------------------------------------	
	vector<double> perrun;          //record vMSD of all atoms in each run for each type (time average)
	vector<double> allrun;          //record vMSD of all runs for each type (atom average)
	vector<double> vmsdPerType;     //record vMSD for each type (run average)
	char filename[64];
	double averagePerRun;
	double averagePerType;
	double totalPerType;

	for(int atomType = 0; atomType < thisNumPerType.getsize(); ++ atomType)
	{
		for(int runnumber = 1; runnumber <= 100; ++ runnumber)
		{
			sprintf(filename, "%d/XDATCAR", runnumber);
			readCoordination(filename);
			modifyCoordination();
			averagePerRun = vibrationalmsd(perrun, atomType, thisNumPerType);
			allrun.push_back(averagePerRun);

			perrun.clear();
		}
		for(int runnumber = 0; runnumber < allrun.size(); ++ runnumber)
		{
			totalPerType += allrun[runnumber];
		}
		averagePerType = totalPerType /allrun.size();
		vmsdPerType.push_back(averagePerType);

		allrun.clear();
	}	

	fp = fopen("vibrationalMSD_information", "w");
	for(int i = 0; i < thisNumPerType.getsize(); ++ i)
	{
		fprintf(fp, "%s\t%.8lf\n", namelist[i].c_str(), vmsdPerType[i]);
	}
	fclose(fp);

	return 0;
}
