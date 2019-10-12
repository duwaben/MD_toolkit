/*
	Name: Pair Correlation Function calculation program 3.0
	Time: 20/06/2018
	Input file: XDATCAR
	Output file: PDF11_8x, PDF22_8x, PDF33_8x, PDF12_8x, PDF13_8x, PDF23_8x.
	Improvement: improve the efficiency 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXBIN 1000 //the half of box length(multiply length by 2 in three dimensions) 
#define delta_r 0.01
#define gap 1
#define Pi 3.141593
#define NOMIT 4000
#define NSTEP 8000
#define NATOM 200

FILE *fp;
double u1[NATOM], v1[NATOM], w1[NATOM];
double LATTICE;
int NATOM_1, NATOM_2, NATOM_3;//three types

double gr[MAXBIN];

//functional function
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

void r(double uc1, double vc1, double wc1, double uc2, double vc2, double wc2, double *rpbc)
{
	double x, y, z;
	double LATTICE_8x = 2 * LATTICE;
	x = uc1 - uc2;
	y = vc1 - vc2;
	z = wc1 - wc2;
	if(x > LATTICE_8x /2) x -= LATTICE_8x; if(x < -LATTICE_8x /2) x += LATTICE_8x;
	if(y > LATTICE_8x /2) y -= LATTICE_8x; if(y < -LATTICE_8x /2) y += LATTICE_8x;
	if(z > LATTICE_8x /2) z -= LATTICE_8x; if(z < -LATTICE_8x /2) z += LATTICE_8x;
	*rpbc = sqrt(x * x + y * y + z * z);
}

void PDF_XX_8x(int numstart1, int numend1)
{
	int i, j;
  	int BIN;
  	int num1 = numend1 - numstart1;
  	double u_8x[num1 * 8], v_8x[num1 * 8], w_8x[num1 * 8];
  	double delr = delta_r;
  	double rpbc;
  	int hist[MAXBIN];
  	double LATTICE_8x = 2 * LATTICE;

	for(j = numstart1; j < numend1; j ++){
		u_8x[j-numstart1] = u1[j];
		v_8x[j-numstart1] = v1[j];
		w_8x[j-numstart1] = w1[j];

		u_8x[j-numstart1+num1] = u1[j] + LATTICE;
		v_8x[j-numstart1+num1] = v1[j];
		w_8x[j-numstart1+num1] = w1[j];

		u_8x[j-numstart1+num1*2] = u1[j];
		v_8x[j-numstart1+num1*2] = v1[j] + LATTICE;
		w_8x[j-numstart1+num1*2] = w1[j];

		u_8x[j-numstart1+num1*3] = u1[j]; 
		v_8x[j-numstart1+num1*3] = v1[j]; 
		w_8x[j-numstart1+num1*3] = w1[j] + LATTICE;

		u_8x[j-numstart1+num1*4] = u1[j] + LATTICE; 
		v_8x[j-numstart1+num1*4] = v1[j] + LATTICE; 
		w_8x[j-numstart1+num1*4] = w1[j]; 

		u_8x[j-numstart1+num1*5] = u1[j] + LATTICE; 
		v_8x[j-numstart1+num1*5] = v1[j]; 
		w_8x[j-numstart1+num1*5] = w1[j] + LATTICE;

		u_8x[j-numstart1+num1*6] = u1[j]; 
		v_8x[j-numstart1+num1*6] = v1[j] + LATTICE; 
		w_8x[j-numstart1+num1*6] = w1[j] + LATTICE; 

		u_8x[j-numstart1+num1*7] = u1[j] + LATTICE; 
		v_8x[j-numstart1+num1*7] = v1[j] + LATTICE; 
		w_8x[j-numstart1+num1*7] = w1[j] + LATTICE; 
	}

  	for(BIN = 0; BIN < MAXBIN; BIN ++)
    	hist[BIN] = 0;
  	for(i = 0; i < num1 * 8; ++ i){
  		for(j = 0; j < num1 * 8; ++ j)
  			if(i != j){
				r(u_8x[i], v_8x[i], w_8x[i], u_8x[j], v_8x[j], w_8x[j], &rpbc);  
  				BIN = floor(rpbc /delr); 
  				if(BIN < MAXBIN){
    				hist[BIN] ++; 
  				}
    		}
  	}

  	for(BIN = 0; BIN < MAXBIN; BIN ++){
   		gr[BIN] = (LATTICE_8x * LATTICE_8x * LATTICE_8x /(num1 * 8) /(num1 * 8)) * (double)hist[BIN] /(4.0 * Pi * (((double)BIN+0.5) * delr) * (((double)BIN+0.5) * delr) * delr);
  	}
}

int main()
{
	fp = fopen("XDATCAR", "r");
	
	void Findword(const char *word);
	
	int i, j, k;
	
	const char *word1 = "Direct configuration";
	char line[512];

	double *u[NATOM], *v[NATOM], *w[NATOM];	
	for(i = 0; i < NATOM; i ++){
		u[i] = (double *)calloc(NSTEP, sizeof(double));
		v[i] = (double *)calloc(NSTEP, sizeof(double));
		w[i] = (double *)calloc(NSTEP, sizeof(double));
	}
	
	//read useful informations form files
	//**********************************************************
	//**********************************************************
	fgets(line, 512, fp); fgets(line, 512, fp);
	fscanf(fp, "%lf", &LATTICE);

	//skip the steps
	for(i = 0; i < NOMIT; i ++){
		Findword(word1);
	}	
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
	fclose(fp);
	//*********************************************************
	//*********************************************************
	

	//*********************************************************
	//*********************************************************
	int step = 0; double sumgr[MAXBIN];
	for(i = 0; i < MAXBIN; i ++)
		sumgr[i] = 0.0;
		
	for(i = 0; i < NSTEP; i = i+gap){	
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];	
		}
		step ++;
		PDF_XX_8x(0, NATOM);
		for(j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j]; 
	}
	//step average
	for(j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step; 
		
	if((fp = fopen("pairCorrelationFunctions/PDFtotal_8x","w")) == NULL){
		printf("Cannot create file PDFtotal."); exit(0);
	}

	for(k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * (k+0.5), sumgr[k]);
	}
	fclose(fp);

	return 0;
}
