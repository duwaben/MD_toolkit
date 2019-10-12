#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>
#include <string>
#include <sys/stat.h>

using namespace std;

#define delta_r 0.01
#define gap 10
#define Pi 3.141593
#define NOMIT 27000
#define NSTEP 3000
#define NATOM 200 //NATOM should not be set in my vision.

#define MAX_L   20      //default  20
#define MAXBIN_S   1201 //MAX_L*MAX_L*3+1

FILE *fp;
double LATTICE;
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];
double u1[NATOM], v1[NATOM], w1[NATOM];
double sumcos[MAXBIN_S];
double histcos[MAXBIN_S];

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

void structureFactors_XX(int numstart1, int numend1)
{
	int BIN;
	double qr;
	double du, dv, dw;

	int countBIN[MAXBIN_S];
	for(int i = 0; i < MAXBIN_S; i ++)
	{
		countBIN[i] = 0;
	}

	int Lu[2*MAX_L+1], Lv[2*MAX_L+1], Lw[2*MAX_L+1];
	Lu[0] = Lv[0] = Lw[0] = 0;
	for(int i = 1; i <= MAX_L; ++ i)
	{
		Lu[i] = Lv[i] = Lw[i] = i;
		Lu[i+MAX_L] = Lv[i+MAX_L] = Lw[i+MAX_L] = -i;
	}

	for(int i = 0; i < MAXBIN_S; i ++)
	{
		sumcos[i] = 0.0;
	}

	for(int i = 0; i < 2*MAX_L+1; i ++)
	{
		for(int j = 0; j < 2*MAX_L+1; j ++)
		{
			for (int k = 0; k < 2*MAX_L+1; k ++)
			{
				BIN = Lu[i]*Lu[i] + Lv[j]*Lv[j] + Lw[k]*Lw[k];
				countBIN[BIN] ++;
				for(int ni = numstart1; ni < numend1; ni ++)
				{
					for(int nj = numstart1; nj < numend1; nj ++)
					{
						du = u1[ni]-u1[nj];
						dv = v1[ni]-v1[nj];
						dw = w1[ni]-w1[nj];

						if(du > LATTICE/2) du -= LATTICE;
						if(du < -LATTICE/2) du += LATTICE;
						if(dv > LATTICE/2) dv -= LATTICE;
						if(dv < -LATTICE/2) dv += LATTICE;
						if(dw > LATTICE/2) dw -= LATTICE;
						if(dw < -LATTICE/2) dw += LATTICE;

						qr = du*Lu[i] + dv*Lv[j] + dw*Lw[k];
						sumcos[BIN] += cos(qr * 2 * Pi /LATTICE);
					}
				}
			}
		}
	}

	for(int i = 0; i < MAXBIN_S; i ++)
	{
		if(countBIN[i] != 0)
		{
			sumcos[i] = sumcos[i] /countBIN[i] /(numend1 - numstart1);
		}
	}
}

void structureFactors_XY(int numstart1, int numend1, int numstart2, int numend2)
{
	int BIN;
	double qr;
	double du, dv, dw;

	int countBIN[MAXBIN_S];
	for(int i = 0; i < MAXBIN_S; i ++)
	{
		countBIN[i] = 0;
	}

	int Lu[2*MAX_L+1], Lv[2*MAX_L+1], Lw[2*MAX_L+1];
	Lu[0] = Lv[0] = Lw[0] = 0;
	for(int i = 1; i <= MAX_L; ++ i)
	{
		Lu[i] = Lv[i] = Lw[i] = i;
		Lu[i+MAX_L] = Lv[i+MAX_L] = Lw[i+MAX_L] = -i;
	}

	for(int i = 0; i < MAXBIN_S; i ++)
	{
		sumcos[i] = 0.0;
	}

	for(int i = 0; i < 2*MAX_L+1; i ++)
	{
		for(int j = 0; j < 2*MAX_L+1; j ++)
		{
			for (int k = 0; k < 2*MAX_L+1; k ++)
			{
				BIN = Lu[i]*Lu[i] + Lv[j]*Lv[j] + Lw[k]*Lw[k];
				countBIN[BIN] ++;
				for(int ni = numstart1; ni < numend1; ni ++)
				{
					for(int nj = numstart2; nj < numend2; nj ++)
					{
						du = u1[ni]-u1[nj];
						dv = v1[ni]-v1[nj];
						dw = w1[ni]-w1[nj];

						if(du > LATTICE/2) du -= LATTICE;
						if(du < -LATTICE/2) du += LATTICE;
						if(dv > LATTICE/2) dv -= LATTICE;
						if(dv < -LATTICE/2) dv += LATTICE;
						if(dw > LATTICE/2) dw -= LATTICE;
						if(dw < -LATTICE/2) dw += LATTICE;

						qr = du*Lu[i] + dv*Lv[j] + dw*Lw[k];
						sumcos[BIN] += cos(qr * 2 * Pi /LATTICE);
					}
				}
			}
		}
	}

	for(int i = 0; i < MAXBIN_S; i ++)
	{
		if(countBIN[i] != 0)
		{
			sumcos[i] = sumcos[i] /countBIN[i] /sqrt((numend1 - numstart1)*(numend2 - numstart2));
		}
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

int main(int argc, char const *argv[])
{
	//namelist
	vector<string> namelist;
	readName(namelist);

	//
	NumPerType thisNumPerType;
	thisNumPerType = readCoordinate("XDATCAR");
	mkdir("structureFactors", 00777);	
	
	int step;
	char filename[128];
	for(int i = 0; i < thisNumPerType.getsize(); i ++)
	{
		for(int j = i; j < thisNumPerType.getsize(); j ++)
		{
			step = 0;
			for(int k = 0; k < MAXBIN_S; k ++)
			{
				sumcos[k] = 0.0;
				histcos[k] = 0.0;
			}

			for(int ni = 0; ni < NSTEP; ni += gap)
			{
				step ++;
				for(int nj = 0; nj < NATOM; nj ++)
				{
					u1[nj] = u[nj][ni];
					v1[nj] = v[nj][ni];
					w1[nj] = w[nj][ni];
				}
				if(i == j)
				{
					structureFactors_XX(thisNumPerType.getmin(i), thisNumPerType.getmax(i));
				}
				else
				{
					structureFactors_XY(thisNumPerType.getmin(i), thisNumPerType.getmax(i), thisNumPerType.getmin(j), thisNumPerType.getmax(j));
				}

				for(int l = 0; l < MAXBIN_S; l ++)
				{
					histcos[l] += sumcos[l];
				}
			}

			for(int k = 0; k < MAXBIN_S; k ++)
			{
				histcos[k] = histcos[k] /step;
			}

			sprintf(filename, "structureFactors/SF_%d%d_%s%s", i+1, j+1, namelist[i].c_str(), namelist[j].c_str());
			fp = fopen(filename, "w");
			fprintf(fp, "X\t%s%s\n", namelist[i].c_str(), namelist[j].c_str());
			for(int k = 1; k < MAXBIN_S; k ++)
			{
				if(histcos[k] != 0)
				{
					fprintf(fp, "%lf\t%lf\n", sqrt(k) * 2 * Pi /LATTICE, histcos[k]);					
				}
			}
			fclose(fp);
		}
	}

	return 0;
}
