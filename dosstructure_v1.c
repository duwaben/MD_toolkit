/*
	Name: POSCAR generator(for dos calculation)_v1
	Author: Frank Chen
	Date: 30/03/2018
	Description: the input file is XDATCAR here and it generate 10 configurations for following dos calculation 
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define NATOM 200
#define NOMIT 4000
#define NSTEP 8000

FILE *fp;

void Findword(char *word){
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

void main()
{  
		
	int i, j;
	char *word1 = "Direct configuration";
	char line[512];
	char line1[512], line6[512], line7[512];
	char string[50];
	double LATTICE;
	int p;
	
	double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];

	fp = fopen("XDATCAR", "r");

	fgets(line1, 512, fp); fgets(line, 512, fp);
	fscanf(fp, "%lf", &LATTICE);
	fgets(line, 512, fp); fgets(line, 512, fp); fgets(line, 512, fp);
	fgets(line6, 512, fp); fgets(line7, 512, fp);	
	
	for(i = 0; i < NOMIT; i ++){
		Findword(word1);
		printf("Read step %d.\n", i+1);
	}	
	for(i = 0; i < NSTEP; i ++){
		Findword(word1);
		fgets(line, 512, fp);
		printf("Read step %d.\n", i+2001);
		for(j = 0; j < NATOM; j ++){
			fscanf(fp, "%lf%lf%lf", &u[j][i], &v[j][i], &w[j][i]); 
		}
	}
	fclose(fp);

	for(p = 1; p <= 10; p ++){
		i = 4000 + p * 400 - 1;
		sprintf(string, "POSCAR%d", p);//good strategy
		FILE *fp = fopen(string, "w");
		fputs(line1, fp); fprintf(fp, "%lf\n1.000000  0.000000  0.000000\n0.000000  1.000000  0.000000\n0.000000  0.000000  1.000000\n", LATTICE);
		fputs(line6, fp); fputs(line7, fp); fprintf(fp, "Direct\n");
		for(j = 0; j < NATOM; j ++){	
			fprintf(fp, "\t%lf\t%lf\t%lf\n", u[j][i], v[j][i], w[j][i]);
		}
		fclose(fp);
	}
}
	

			