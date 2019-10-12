#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;

#define NATOM 200       //total atom number
#define NOMIT 4000 		//the number of omitted steps
#define NSTEP 8000		//the number of counting steps
#define gap 1			//the interval of counting
#define Pi 3.141592653

FILE *fp;
double LATTICE;
int NATOM_1, NATOM_2, NATOM_3; //three types
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];	
double u1[NATOM], v1[NATOM], w1[NATOM];

void Findword(const char *word)
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

void readCoordinate(const char *name)
{
	fp = fopen(name, "r");
	
	int i, j;
	
	const char *word1 = "Direct configuration";
	char line[512];

	//read useful informations form files
	//**********************************************************
	//**********************************************************
	fgets(line, 512, fp); fgets(line, 512, fp);
	fscanf(fp, "%lf", &LATTICE);
	
	for(i = 0; i < 4; i ++)fgets(line, 512, fp);	
	fscanf(fp, "%d%d%d", &NATOM_1, &NATOM_2, &NATOM_3);
	
	//omit the steps
	for(i = 0; i < NOMIT; i ++){
		Findword(word1);
	}
	//read the steps	
	for(i = 0; i < NSTEP; i ++){
		Findword(word1);
		fgets(line, 512, fp);
		for(j = 0; j < NATOM; j ++){
			//read the last 6000 steps
			fscanf(fp, "%lf%lf%lf", &u[j][i], &v[j][i], &w[j][i]);
			u[j][i] *= LATTICE;
			v[j][i] *= LATTICE;
			w[j][i] *= LATTICE;
		} 
	}
	printf("Coordinate-reading part is over.\n");
	//read the number of each atom type and the lattice constant
	printf("LATTICE=%lf\tNATOM_1=%d\tNATOM_2=%d\tNATOM_3=%d\n", LATTICE, NATOM_1, NATOM_2, NATOM_3);
	fclose(fp);
	//*********************************************************
	//*********************************************************
}

int main(int argc, char const *argv[])
{
	double CUTOFF;
	CUTOFF = atof(argv[1]);

	readCoordinate("XDATCAR");

	int temp;
	int hist[16];
	for(int i = 0; i < 16; ++i)
		hist[i] = 0;

	for(int i = 0; i < NSTEP; ++ i)
	{
		for(int j = 0; j < NATOM; ++ j)
		{
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}

		for(int j = 0; j < NATOM_1; ++ j)//Ni atoms
		{
			temp = 0;
			for(int k = NATOM_1; k < NATOM_1+NATOM_2; ++ k)
			{
				if(r(u1[j], v1[j], w1[j], u1[k], v1[k], w1[k]) < CUTOFF)
					temp ++;
			}
			hist[temp] ++;
		}
	}

	int totalhist = 0;
	for(int i = 0; i < 16; ++ i)
	{
		totalhist += hist[i];
	}

	double histratio[16];
	for(int i = 0; i < 16; ++ i)
	{
		histratio[i] = (double)hist[i]/(double)totalhist;
	}

	fp = fopen("pairCorrelationFunctions/CNdistribution_NiAl", "w");
	for(int i = 0; i < 16; ++i)
	{
		fprintf(fp, "%d\t%lf\n", i, histratio[i]);
	}
	fclose(fp);

	return 0;
}

