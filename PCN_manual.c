/*
	Name: Pair Coordination Number_v1
	Time: 02/04/2018
	Input file: XDATCAR, CUTOFF
	Output file: PDF11_8x, PDF22_8x, PDF33_8x, PDF12_8x, PDF13_8x, PDF23_8x.
	Improvement: calculate the partial coordination number
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#define MAXBIN 1000 //the half of box length(multiply length by 2 in three dimensions)
#define delta_r 0.01
#define gap 1
#define Pi 3.141593
#define NOMIT 4000 //temperature is so high
#define NSTEP 8000
#define NATOM 200

FILE *fp;
double LATTICE;
int NATOM_1, NATOM_2, NATOM_3;
double CUTOFF11, CUTOFF12, CUTOFF13, CUTOFF22, CUTOFF23, CUTOFF33;
double PEAK11, PEAK12, PEAK13, PEAK22, PEAK23, PEAK33;
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];
double u1[NATOM], v1[NATOM], w1[NATOM];

double gr[MAXBIN];

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

void PDF_XY_8x(int numstart1, int numend1, int numstart2, int numend2)
{
  	int i, j;
  	int BIN;
  	int num1 = numend1 - numstart1;
  	int num2 = numend2 - numstart2;
 	double u_8x[num1*8+num2*8], v_8x[num1*8+num2*8], w_8x[num1*8+num2*8];
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

	for(j = numstart2; j < numend2; j ++){
		u_8x[j-numstart2+num1*8] = u1[j];
		v_8x[j-numstart2+num1*8] = v1[j];
		w_8x[j-numstart2+num1*8] = w1[j];

		u_8x[j-numstart2+num1*8+num2] = u1[j] + LATTICE;
		v_8x[j-numstart2+num1*8+num2] = v1[j];
		w_8x[j-numstart2+num1*8+num2] = w1[j];

		u_8x[j-numstart2+num1*8+num2*2] = u1[j];
		v_8x[j-numstart2+num1*8+num2*2] = v1[j] + LATTICE;
		w_8x[j-numstart2+num1*8+num2*2] = w1[j];

		u_8x[j-numstart2+num1*8+num2*3] = u1[j];
		v_8x[j-numstart2+num1*8+num2*3] = v1[j];
		w_8x[j-numstart2+num1*8+num2*3] = w1[j] + LATTICE;

		u_8x[j-numstart2+num1*8+num2*4] = u1[j] + LATTICE;
		v_8x[j-numstart2+num1*8+num2*4] = v1[j] + LATTICE;
		w_8x[j-numstart2+num1*8+num2*4] = w1[j];

		u_8x[j-numstart2+num1*8+num2*5] = u1[j] + LATTICE;
		v_8x[j-numstart2+num1*8+num2*5] = v1[j];
		w_8x[j-numstart2+num1*8+num2*5] = w1[j] + LATTICE;

		u_8x[j-numstart2+num1*8+num2*6] = u1[j];
		v_8x[j-numstart2+num1*8+num2*6] = v1[j] + LATTICE;
		w_8x[j-numstart2+num1*8+num2*6] = w1[j] + LATTICE;

		u_8x[j-numstart2+num1*8+num2*7] = u1[j] + LATTICE;
		v_8x[j-numstart2+num1*8+num2*7] = v1[j] + LATTICE;
		w_8x[j-numstart2+num1*8+num2*7] = w1[j] + LATTICE;
	}

  	for(BIN = 0; BIN < MAXBIN; BIN ++)
    	hist[BIN] = 0;
  	for(i = 0; i < num1*8; ++ i){
  		for(j = num1*8; j < num1*8+num2*8; ++ j){
			r(u_8x[i], v_8x[i], w_8x[i], u_8x[j], v_8x[j], w_8x[j], &rpbc);
  			BIN = floor(rpbc /delr);
  			if(BIN < MAXBIN)
    			hist[BIN] ++;
    	}
  	}
  	for(BIN = 0; BIN < MAXBIN; BIN ++){
   		gr[BIN] = (LATTICE_8x * LATTICE_8x * LATTICE_8x /(num1 * 8) /(num2 * 8)) * (double)hist[BIN] /(4.0 * Pi * (((double)BIN+0.5) * delr) * (((double)BIN+0.5) * delr) * delr);
  	}
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

	//skip the steps
	for(i = 0; i < NOMIT; i ++){
		findWord(word1);
	}
	for(i = 0; i < NSTEP; i ++){
		findWord(word1);
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

void readCutoff(const char *name)
{
	fp = fopen(name, "r");

	const char *cutoffword1 = "CUTOFF11:";
	const char *cutoffword2 = "CUTOFF12:";
	const char *cutoffword3 = "CUTOFF13:";
	const char *cutoffword4 = "CUTOFF22:";
	const char *cutoffword5 = "CUTOFF23:";
	const char *cutoffword6 = "CUTOFF33:";
	const char *peakword1 = "PEAK11:";
	const char *peakword2 = "PEAK12:";
	const char *peakword3 = "PEAK13:";
	const char *peakword4 = "PEAK22:";
	const char *peakword5 = "PEAK23:";
	const char *peakword6 = "PEAK33:";

	findWord(cutoffword1);
	fscanf(fp, "%lf", &CUTOFF11);
	findWord(cutoffword2);
	fscanf(fp, "%lf", &CUTOFF12);
	findWord(cutoffword3);
	fscanf(fp, "%lf", &CUTOFF13);
	findWord(cutoffword4);
	fscanf(fp, "%lf", &CUTOFF22);
	findWord(cutoffword5);
	fscanf(fp, "%lf", &CUTOFF23);
	findWord(cutoffword6);
	fscanf(fp, "%lf", &CUTOFF33);

	findWord(peakword1);
	fscanf(fp, "%lf", &PEAK11);
	findWord(peakword2);
	fscanf(fp, "%lf", &PEAK12);
	findWord(peakword3);
	fscanf(fp, "%lf", &PEAK13);
	findWord(peakword4);
	fscanf(fp, "%lf", &PEAK22);
	findWord(peakword5);
	fscanf(fp, "%lf", &PEAK23);
	findWord(peakword6);
	fscanf(fp, "%lf", &PEAK33);

	fclose(fp);
}

int main()
{
	int i, j, k;
	double pCN_11, pCN_12, pCN_13, pCN_21, pCN_22, pCN_23, pCN_31, pCN_32, pCN_33;
	readCoordinate("XDATCAR");
	readCutoff("pairCorrelationFunctions/CUTOFFPEAK_information");

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
		PDF_XX_8x(0, NATOM_1);
		for(j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j];
	}
	//step average
	for(j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step;

	//pCN
	double sumpCN = 0;
	int CUTOFFBIN = floor(CUTOFF11/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_11 = sumpCN * 4 * Pi * NATOM_1 /(LATTICE * LATTICE * LATTICE);
    //----------------------------------------------------------

	step = 0;
	for(i = 0; i < MAXBIN; i ++)
		sumgr[i] = 0.0;

	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		step ++;
		PDF_XX_8x(NATOM_1, NATOM_1+NATOM_2);
		for(j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j];
	}
	//step average
	for(j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step;

	//pCN
	sumpCN = 0;
	CUTOFFBIN = floor(CUTOFF22/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_22 = sumpCN * 4 * Pi * NATOM_2 /(LATTICE * LATTICE * LATTICE);

    //---------------------------------------------------------

	step = 0;
	for(i = 0; i < MAXBIN; i ++)
		sumgr[i] = 0.0;

	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		step ++;
		PDF_XX_8x(NATOM_1+NATOM_2, NATOM);
		for(j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j];
	}
	//step average
	for(j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step;

	//pCN
	sumpCN = 0;
	CUTOFFBIN = floor(CUTOFF33/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_33 = sumpCN * 4 * Pi * NATOM_3 /(LATTICE * LATTICE * LATTICE);
	//--------------------------------------------------------

	//PDF_XY
	step = 0;
	for(i = 0; i < MAXBIN; i ++)
		sumgr[i] = 0.0;

	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		step ++;
		PDF_XY_8x(0, NATOM_1, NATOM_1, NATOM_1+NATOM_2);
		for(j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j];
	}
	//step average
	for(j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step;

	//pCN
	sumpCN = 0;
	CUTOFFBIN = floor(CUTOFF12/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_12 = sumpCN * 4 * Pi * NATOM_2 /(LATTICE * LATTICE * LATTICE);
	pCN_21 = sumpCN * 4 * Pi * NATOM_1 /(LATTICE * LATTICE * LATTICE);

	//----------------------------------------------------------


	step = 0;
	for(i = 0; i < MAXBIN; i ++)
		sumgr[i] = 0.0;

	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		step ++;
		PDF_XY_8x(0, NATOM_1, NATOM_1+NATOM_2, NATOM);
		for(j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j];
	}
	//step average
	for(j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step;

	//pCN
	sumpCN = 0;
	CUTOFFBIN = floor(CUTOFF13/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_13 = sumpCN * 4 * Pi * NATOM_3 /(LATTICE * LATTICE * LATTICE);
	pCN_31 = sumpCN * 4 * Pi * NATOM_1 /(LATTICE * LATTICE * LATTICE);

	//----------------------------------------------------------


	step = 0;
	for(i = 0; i < MAXBIN; i ++)
		sumgr[i] = 0.0;

	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		step ++;
		PDF_XY_8x(NATOM_1, NATOM_1+NATOM_2, NATOM_1+NATOM_2, NATOM);
		for(j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j];
	}
	//step average
	for(j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step;

	//pCN
	sumpCN = 0;
	CUTOFFBIN = floor(CUTOFF23/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_23 = sumpCN * 4 * Pi * NATOM_3 /(LATTICE * LATTICE * LATTICE);
	pCN_32 = sumpCN * 4 * Pi * NATOM_2 /(LATTICE * LATTICE * LATTICE);


	//***********************************************************************
	//***********************************************************************
	//record the information of pCN.
	fp = fopen("pairCorrelationFunctions/pCN_information", "w");

	fprintf(fp, "pCN_11: %lf\n", pCN_11);
	fprintf(fp, "pCN_12: %lf\n", pCN_12);
	fprintf(fp, "pCN_13: %lf\n", pCN_13);
	fprintf(fp, "pCN_21: %lf\n", pCN_21);
	fprintf(fp, "pCN_22: %lf\n", pCN_22);
	fprintf(fp, "pCN_23: %lf\n", pCN_23);
	fprintf(fp, "pCN_31: %lf\n", pCN_31);
	fprintf(fp, "pCN_32: %lf\n", pCN_32);
	fprintf(fp, "pCN_33: %lf\n", pCN_33);

	fclose(fp);
	//***********************************************************************
	//***********************************************************************

	//***********************************************************************
	//***********************************************************************
	//record the information of Warren-Cowley parameter.
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	double Z1, Z2, Z3;
	fp = fopen("pairCorrelationFunctions/WC_information", "w");
	Z1 = pCN_11 + pCN_12 + pCN_13;
	Z2 = pCN_21 + pCN_22 + pCN_23;
	Z3 = pCN_31 + pCN_32 + pCN_33;
	if(NATOM_3 == 0){
		Z1 = pCN_11 + pCN_12;
		Z2 = pCN_21 + pCN_22;		
	}

	a11 = 1.0 - pCN_11/((double)NATOM_1/NATOM * Z1);
	a12 = 1.0 - pCN_12/((double)NATOM_2/NATOM * Z1);
	a13 = 1.0 - pCN_13/((double)NATOM_3/NATOM * Z1);
	a21 = 1.0 - pCN_21/((double)NATOM_1/NATOM * Z2);
	a22 = 1.0 - pCN_22/((double)NATOM_2/NATOM * Z2);
	a23 = 1.0 - pCN_23/((double)NATOM_3/NATOM * Z2);
	a31 = 1.0 - pCN_31/((double)NATOM_1/NATOM * Z3);
	a32 = 1.0 - pCN_32/((double)NATOM_2/NATOM * Z3);
	a33 = 1.0 - pCN_33/((double)NATOM_3/NATOM * Z3);

	fprintf(fp, "a11: %lf\n", a11);
	fprintf(fp, "a12: %lf\n", a12);
	fprintf(fp, "a13: %lf\n", a13);
	fprintf(fp, "a21: %lf\n", a21);
	fprintf(fp, "a22: %lf\n", a22);
	fprintf(fp, "a23: %lf\n", a23);
	fprintf(fp, "a31: %lf\n", a31);
	fprintf(fp, "a32: %lf\n", a32);
	fprintf(fp, "a33: %lf\n", a33);

	fclose(fp);
	//***********************************************************************
	//***********************************************************************
	
	//***********************************************************************
	//***********************************************************************
	//record the compactness
	double sigma11, sigma12, sigma13, sigma21, sigma22, sigma23, sigma31, sigma32, sigma33; 
	
	sigma11 = pCN_11/(4 * Pi * PEAK11 * PEAK11);
	sigma12 = pCN_12/(4 * Pi * PEAK12 * PEAK12);
	sigma13 = pCN_13/(4 * Pi * PEAK13 * PEAK13);
	sigma21 = pCN_21/(4 * Pi * PEAK12 * PEAK12);
	sigma22 = pCN_22/(4 * Pi * PEAK22 * PEAK22);
	sigma23 = pCN_23/(4 * Pi * PEAK23 * PEAK23);
	sigma31 = pCN_31/(4 * Pi * PEAK13 * PEAK13);
	sigma32 = pCN_32/(4 * Pi * PEAK23 * PEAK23);
	sigma33 = pCN_33/(4 * Pi * PEAK33 * PEAK33);

	fp = fopen("pairCorrelationFunctions/COMPACTNESS_information", "w");

	fprintf(fp, "sigma11: %lf\n", sigma11);
	fprintf(fp, "sigma12: %lf\n", sigma12);
	fprintf(fp, "sigma13: %lf\n", sigma13);
	fprintf(fp, "sigma21: %lf\n", sigma21);
	fprintf(fp, "sigma22: %lf\n", sigma22);
	fprintf(fp, "sigma23: %lf\n", sigma23);
	fprintf(fp, "sigma31: %lf\n", sigma31);
	fprintf(fp, "sigma32: %lf\n", sigma32);
	fprintf(fp, "sigma33: %lf\n", sigma33);

	fclose(fp);
	//***********************************************************************
	//***********************************************************************

	return 0;
}
