/*
	Name: Common Neighbour Analysis 1.0
	Date: 01/03/2018
	Input file: XDATCAR, CUTOFF
	Output file: CNA_information.
	Improvement:
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define gap 1
#define Pi 3.141593
#define NOMIT 4000
#define NSTEP 8000
#define NATOM 200

FILE *fp;
double LATTICE;
int NATOM_1, NATOM_2, NATOM_3;
double CUTOFF11, CUTOFF12, CUTOFF13, CUTOFF22, CUTOFF23, CUTOFF33;
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];
double u1[NATOM], v1[NATOM], w1[NATOM];

int hist1661, hist1551, hist1541, hist1441, hist1431, hist142x, histbonds;

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
	x = uc1 - uc2;
	y = vc1 - vc2;
	z = wc1 - wc2;
	if(x > LATTICE /2) x -= LATTICE; if(x < -LATTICE /2) x += LATTICE;
	if(y > LATTICE /2) y -= LATTICE; if(y < -LATTICE /2) y += LATTICE;
	if(z > LATTICE /2) z -= LATTICE; if(z < -LATTICE /2) z += LATTICE;
	*rpbc = sqrt(x * x + y * y + z * z);
}

void commonNeighbourAnalysis_XX(int startnum1, int endnum1, double R_cut)
{
	int i, j, ni, nj;
	int neighbourcount;
	int *array;
	double rab;

	hist1661 = hist1551 = hist1541 = hist1441 = hist1431 = hist142x = 0;
 	histbonds = 0;
 	for(i = startnum1; i < endnum1; i ++){
 		for(j = i+1; j < endnum1; j ++){
 			
			//find out if the two atoms are bonding
			r(u1[i], v1[i], w1[i], u1[j], v1[j], w1[j], &rab);
			if(rab < R_cut)
				histbonds ++;
			else
				continue;

			//********************************************************************
			//********************************************************************
			//figure out the number of neighbours
			neighbourcount = 0;
			array = (int *)malloc(15 * sizeof(int));
			for(ni = 0; ni < NATOM; ni ++){
				if(ni != i && ni != j){
					r(u1[ni], v1[ni], w1[ni], u1[i], v1[i], w1[i], &rab);
					if(rab < R_cut){
						r(u1[ni], v1[ni], w1[ni], u1[j], v1[j], w1[j], &rab);
						if(rab < R_cut){
							array[neighbourcount] = ni;
							neighbourcount ++;
						}
					}	 					
				}
			}
			array = (int *)realloc(array, neighbourcount * sizeof(int));
 				
			//*******************************************************************
			//*******************************************************************
			//figure out the number of bonding neighbours
			int sumbond, productbond;
			int bondspernode;

			if(neighbourcount == 6){
				sumbond = 0;
				productbond = 1;
				
				for(ni = 0; ni < neighbourcount; ni ++){
					bondspernode = 0;
					for(nj = 0; nj < neighbourcount; nj ++)
						if(ni != nj){
	 						r(u1[array[ni]], v1[array[ni]], w1[array[ni]], u1[array[nj]], v1[array[nj]], w1[array[nj]], &rab);
	 						if(rab < R_cut){
	 							sumbond ++;
	 							bondspernode ++;
	 						}
						}
					productbond *= bondspernode;
				}

				if(sumbond == 12 && productbond == 64){hist1661++;}
				if(sumbond == 14 && productbond == 144){hist1661++;}
				if(sumbond == 16 && productbond == 288){hist1661++;}
				if(sumbond == 16 && productbond == 324){hist1661++;}
				if(sumbond == 18 && productbond == 512){hist1661++;}
				if(sumbond == 18 && productbond == 576){hist1661++;}//
				if(sumbond == 18 && productbond == 648){hist1661++;}
				if(sumbond == 20 && productbond == 1152){hist1661++;}//
				if(sumbond == 20 && productbond == 1296){hist1661++;}
				if(sumbond == 22 && productbond == 2304){hist1661++;}
				if(sumbond == 24 && productbond == 4096){hist1661++;}
			}

			//-------------------------------------------------------------	
			if(neighbourcount == 5){
				sumbond = 0;
				productbond = 1;
				
				for(ni = 0; ni < neighbourcount; ni ++){
					bondspernode = 0;
					for(nj = 0; nj < neighbourcount; nj ++)
						if(ni != nj){
	 						r(u1[array[ni]], v1[array[ni]], w1[array[ni]], u1[array[nj]], v1[array[nj]], w1[array[nj]], &rab);
	 						if(rab < R_cut){
	 							sumbond ++;
	 							bondspernode ++;
 							}
						}
					productbond *= bondspernode;
				}

				if(sumbond == 10 && productbond == 32){hist1551++;}
				if(sumbond == 12 && productbond == 72){hist1551++;}
				if(sumbond == 14 && productbond == 144){hist1551++;}
				if(sumbond == 14 && productbond == 162){hist1551++;}
				if(sumbond == 16 && productbond == 288){hist1551++;}
				if(sumbond == 16 && productbond == 324){hist1551++;}
				if(sumbond == 18 && productbond == 576){hist1551++;}
				if(sumbond == 20 && productbond == 1024){hist1551++;}

				if(sumbond == 8  && productbond == 8){hist1541++;}
				if(sumbond == 10 && productbond == 18){hist1541++;}
				if(sumbond == 10 && productbond == 24){hist1541++;}
				if(sumbond == 12 && productbond == 48){hist1541++;}//
				if(sumbond == 12 && productbond == 54){hist1541++;}
				if(sumbond == 12 && productbond == 64){hist1541++;}
				if(sumbond == 14 && productbond == 108){hist1541++;}//

			}

			//------------------------------------------------------------------------
			if(neighbourcount == 4){
				sumbond = 0;
				productbond = 1;
				
				for(ni = 0; ni < neighbourcount; ni ++){
					bondspernode = 0;
					for(nj = 0; nj < neighbourcount; nj ++)
						if(ni != nj){
	 						r(u1[array[ni]], v1[array[ni]], w1[array[ni]], u1[array[nj]], v1[array[nj]], w1[array[nj]], &rab);
	 						if(rab < R_cut){
	 							sumbond ++;
	 							bondspernode ++;
	 						}
						}
					productbond *= bondspernode;
				}

				if(sumbond == 8 && productbond == 16){hist1441++;}
				if(sumbond == 10 && productbond == 36){hist1441++;}
				if(sumbond == 12 && productbond == 81){hist1441++;}
				if(sumbond == 6 && productbond == 4){hist1431++;}
				if(sumbond == 8 && productbond == 12){hist1431++;}
				if(sumbond == 4 && productbond == 1){hist142x++;}
				if(sumbond == 4 && productbond == 0){hist142x++;}
				if(sumbond == 6 && productbond == 0){hist142x++;}

			}

			free(array);
 			
 		}
 	}
}

void commonNeighbourAnalysis_XY(int startnum1, int endnum1, int startnum2, int endnum2, double R_cut)
{
	int i, j, ni, nj;
	int neighbourcount;
	int *array;
	double rab;

	hist1661 = hist1551 = hist1541 = hist1441 = hist1431 = hist142x = 0;
 	histbonds = 0;
 	for(i = startnum1; i < endnum1; i ++){
 		for(j = startnum2; j < endnum2; j ++){
 			
 			//find out if the two atoms are bonding
 			r(u1[i], v1[i], w1[i], u1[j], v1[j], w1[j], &rab);
			if(rab < R_cut)
				histbonds ++;
			else
				continue;

			//********************************************************************
			//********************************************************************
			//figure out the number of neighbours
			neighbourcount = 0;
			array = (int *)malloc(15 * sizeof(int));
			for(ni = 0; ni < NATOM; ni ++){
				if(ni != i && ni != j){
					r(u1[ni], v1[ni], w1[ni], u1[i], v1[i], w1[i], &rab);
					if(rab < R_cut){
						r(u1[ni], v1[ni], w1[ni], u1[j], v1[j], w1[j], &rab);
						if(rab < R_cut){
							array[neighbourcount] = ni;
							neighbourcount ++;
						}
					}	 					
				}
			}
			array = (int *)realloc(array, neighbourcount * sizeof(int));
			
			//*******************************************************************
			//*******************************************************************
			//figure out the number of bonding neighbours
			int sumbond, productbond;
			int bondspernode;

			if(neighbourcount == 6){
				sumbond = 0;
				productbond = 1;
					
				for(ni = 0; ni < neighbourcount; ni ++){
					bondspernode = 0;
					for(nj = 0; nj < neighbourcount; nj ++)
						if(ni != nj){
	 					r(u1[array[ni]], v1[array[ni]], w1[array[ni]], u1[array[nj]], v1[array[nj]], w1[array[nj]], &rab);
	 						if(rab < R_cut){
	 							sumbond ++;
	 							bondspernode ++;
	 						}
						}
						productbond *= bondspernode;
					}

				if(sumbond == 12 && productbond == 64){hist1661++;}
				if(sumbond == 14 && productbond == 144){hist1661++;}
				if(sumbond == 16 && productbond == 288){hist1661++;}
				if(sumbond == 16 && productbond == 324){hist1661++;}
				if(sumbond == 18 && productbond == 512){hist1661++;}
				if(sumbond == 18 && productbond == 576){hist1661++;}//
				if(sumbond == 18 && productbond == 648){hist1661++;}
				if(sumbond == 20 && productbond == 1152){hist1661++;}//
				if(sumbond == 20 && productbond == 1296){hist1661++;}
				if(sumbond == 22 && productbond == 2304){hist1661++;}
				if(sumbond == 24 && productbond == 4096){hist1661++;}
			}

			//-------------------------------------------------------------	
 			if(neighbourcount == 5){
				sumbond = 0;
				productbond = 1;
 					
				for(ni = 0; ni < neighbourcount; ni ++){
					bondspernode = 0;
					for(nj = 0; nj < neighbourcount; nj ++)
						if(ni != nj){
 							r(u1[array[ni]], v1[array[ni]], w1[array[ni]], u1[array[nj]], v1[array[nj]], w1[array[nj]], &rab);
	 						if(rab < R_cut){
	 							sumbond ++;
	 							bondspernode ++;
	 						}
						}
					productbond *= bondspernode;
				}

				if(sumbond == 10 && productbond == 32){hist1551++;}
				if(sumbond == 12 && productbond == 72){hist1551++;}
				if(sumbond == 14 && productbond == 144){hist1551++;}
				if(sumbond == 14 && productbond == 162){hist1551++;}
				if(sumbond == 16 && productbond == 288){hist1551++;}
				if(sumbond == 16 && productbond == 324){hist1551++;}
				if(sumbond == 18 && productbond == 576){hist1551++;}
				if(sumbond == 20 && productbond == 1024){hist1551++;}

				if(sumbond == 8  && productbond == 8){hist1541++;}
				if(sumbond == 10 && productbond == 18){hist1541++;}
				if(sumbond == 10 && productbond == 24){hist1541++;}
				if(sumbond == 12 && productbond == 48){hist1541++;}//
				if(sumbond == 12 && productbond == 54){hist1541++;}
				if(sumbond == 12 && productbond == 64){hist1541++;}
				if(sumbond == 14 && productbond == 108){hist1541++;}//
 			}

			//------------------------------------------------------------------------
			if(neighbourcount == 4){
				sumbond = 0;
				productbond = 1;
 					
				for(ni = 0; ni < neighbourcount; ni ++){
					bondspernode = 0;
					for(nj = 0; nj < neighbourcount; nj ++)
						if(ni != nj){
	 						r(u1[array[ni]], v1[array[ni]], w1[array[ni]], u1[array[nj]], v1[array[nj]], w1[array[nj]], &rab);
	 						if(rab < R_cut){
	 							sumbond ++;
	 							bondspernode ++;
	 						}
						}
					productbond *= bondspernode;
				}

				if(sumbond == 8 && productbond == 16){hist1441++;}
				if(sumbond == 10 && productbond == 36){hist1441++;}
				if(sumbond == 12 && productbond == 81){hist1441++;}
				if(sumbond == 6 && productbond == 4){hist1431++;}
				if(sumbond == 8 && productbond == 12){hist1431++;}
				if(sumbond == 4 && productbond == 1){hist142x++;}
				if(sumbond == 4 && productbond == 0){hist142x++;}
				if(sumbond == 6 && productbond == 0){hist142x++;}

 			}

			free(array);
 		}
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

	fclose(fp);
}

int main(){
	int i, j, k;

	readCoordinate("XDATCAR");
	readCutoff("pairCorrelationFunctions/CUTOFFPEAK_information");	

	//****************************************************************
	//****************************************************************
	int sum1661, sum1551, sum1541, sum1441, sum1431, sum142x, sumbonds;
	double frac1661, frac1551, frac1541, frac1441, frac1431, frac142x;

	//-------------------------------------------------------------
	sum1661 = sum1551 = sum1541 = sum1441 = sum1431 = sum142x = sumbonds = 0;
	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		commonNeighbourAnalysis_XX(0, NATOM_1, CUTOFF11);
		sum1661 += hist1661;
		sum1551 += hist1551;
		sum1541 += hist1541;
		sum1441 += hist1441;
		sum1431 += hist1431;
		sum142x += hist142x;
		sumbonds += histbonds;
	}
	frac1661 = (double)sum1661/sumbonds;
	frac1551 = (double)sum1551/sumbonds;
	frac1541 = (double)sum1541/sumbonds;
	frac1441 = (double)sum1441/sumbonds;
	frac1431 = (double)sum1431/sumbonds;
	frac142x = (double)sum142x/sumbonds;

	fp = fopen("pairCorrelationFunctions/CNA_information", "w");
	fprintf(fp, "\t\t1661\t1551\t1541\t1441\t1431\t142x\n");
	fprintf(fp, "bond11\t");
	fprintf(fp, "%lf\t", frac1661);
	fprintf(fp, "%lf\t", frac1551);
	fprintf(fp, "%lf\t", frac1541);
	fprintf(fp, "%lf\t", frac1441);
	fprintf(fp, "%lf\t", frac1431);
	fprintf(fp, "%lf\t", frac142x);
	fclose(fp);


	//-----------------------------------------------------
	sum1661 = sum1551 = sum1541 = sum1441 = sum1431 = sum142x = sumbonds = 0;
	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		commonNeighbourAnalysis_XX(NATOM_1, NATOM_1+NATOM_2, CUTOFF22);
		sum1661 += hist1661;
		sum1551 += hist1551;
		sum1541 += hist1541;
		sum1441 += hist1441;
		sum1431 += hist1431;
		sum142x += hist142x;
		sumbonds += histbonds;
	}
	frac1661 = (double)sum1661/sumbonds;
	frac1551 = (double)sum1551/sumbonds;
	frac1541 = (double)sum1541/sumbonds;
	frac1441 = (double)sum1441/sumbonds;
	frac1431 = (double)sum1431/sumbonds;
	frac142x = (double)sum142x/sumbonds;

	fp = fopen("pairCorrelationFunctions/CNA_information", "a+");
	fprintf(fp, "\nbond22\t");
	fprintf(fp, "%lf\t", frac1661);
	fprintf(fp, "%lf\t", frac1551);
	fprintf(fp, "%lf\t", frac1541);
	fprintf(fp, "%lf\t", frac1441);
	fprintf(fp, "%lf\t", frac1431);
	fprintf(fp, "%lf\t", frac142x);
	fclose(fp);

	//------------------------------------------------------------
	sum1661 = sum1551 = sum1541 = sum1441 = sum1431 = sum142x = sumbonds = 0;
	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		commonNeighbourAnalysis_XX(NATOM_1+NATOM_2, NATOM, CUTOFF33);
		sum1661 += hist1661;
		sum1551 += hist1551;
		sum1541 += hist1541;
		sum1441 += hist1441;
		sum1431 += hist1431;
		sum142x += hist142x;
		sumbonds += histbonds;
	}
	frac1661 = (double)sum1661/sumbonds;
	frac1551 = (double)sum1551/sumbonds;
	frac1541 = (double)sum1541/sumbonds;
	frac1441 = (double)sum1441/sumbonds;
	frac1431 = (double)sum1431/sumbonds;
	frac142x = (double)sum142x/sumbonds;

	fp = fopen("pairCorrelationFunctions/CNA_information", "a+");
	fprintf(fp, "\nbond33\t");
	fprintf(fp, "%lf\t", frac1661);
	fprintf(fp, "%lf\t", frac1551);
	fprintf(fp, "%lf\t", frac1541);
	fprintf(fp, "%lf\t", frac1441);
	fprintf(fp, "%lf\t", frac1431);
	fprintf(fp, "%lf\t", frac142x);
	fclose(fp);
	//----------------------------------------------------------------

	sum1661 = sum1551 = sum1541 = sum1441 = sum1431 = sum142x = sumbonds = 0;
	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		commonNeighbourAnalysis_XY(0, NATOM_1, NATOM_1, NATOM_1+NATOM_2, CUTOFF12);
		sum1661 += hist1661;
		sum1551 += hist1551;
		sum1541 += hist1541;
		sum1441 += hist1441;
		sum1431 += hist1431;
		sum142x += hist142x;
		sumbonds += histbonds;
	}
	frac1661 = (double)sum1661/sumbonds;
	frac1551 = (double)sum1551/sumbonds;
	frac1541 = (double)sum1541/sumbonds;
	frac1441 = (double)sum1441/sumbonds;
	frac1431 = (double)sum1431/sumbonds;
	frac142x = (double)sum142x/sumbonds;

	fp = fopen("pairCorrelationFunctions/CNA_information", "a+");
	fprintf(fp, "\nbond12\t");
	fprintf(fp, "%lf\t", frac1661);
	fprintf(fp, "%lf\t", frac1551);
	fprintf(fp, "%lf\t", frac1541);
	fprintf(fp, "%lf\t", frac1441);
	fprintf(fp, "%lf\t", frac1431);
	fprintf(fp, "%lf\t", frac142x);
	fclose(fp);
	//----------------------------------------------------------------

	sum1661 = sum1551 = sum1541 = sum1441 = sum1431 = sum142x = sumbonds = 0;
	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		commonNeighbourAnalysis_XY(0, NATOM_1, NATOM_1+NATOM_2, NATOM, CUTOFF13);
		sum1661 += hist1661;
		sum1551 += hist1551;
		sum1541 += hist1541;
		sum1441 += hist1441;
		sum1431 += hist1431;
		sum142x += hist142x;
		sumbonds += histbonds;
	}
	frac1661 = (double)sum1661/sumbonds;
	frac1551 = (double)sum1551/sumbonds;
	frac1541 = (double)sum1541/sumbonds;
	frac1441 = (double)sum1441/sumbonds;
	frac1431 = (double)sum1431/sumbonds;
	frac142x = (double)sum142x/sumbonds;

	fp = fopen("pairCorrelationFunctions/CNA_information", "a+");
	fprintf(fp, "\nbond13\t");
	fprintf(fp, "%lf\t", frac1661);
	fprintf(fp, "%lf\t", frac1551);
	fprintf(fp, "%lf\t", frac1541);
	fprintf(fp, "%lf\t", frac1441);
	fprintf(fp, "%lf\t", frac1431);
	fprintf(fp, "%lf\t", frac142x);
	fclose(fp);
	//----------------------------------------------------------------

	sum1661 = sum1551 = sum1541 = sum1441 = sum1431 = sum142x = sumbonds = 0;
	for(i = 0; i < NSTEP; i = i+gap){
		for(j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		commonNeighbourAnalysis_XY(NATOM_1, NATOM_1+NATOM_2, NATOM_1+NATOM_2, NATOM, CUTOFF23);
		sum1661 += hist1661;
		sum1551 += hist1551;
		sum1541 += hist1541;
		sum1441 += hist1441;
		sum1431 += hist1431;
		sum142x += hist142x;
		sumbonds += histbonds;
	}
	frac1661 = (double)sum1661/sumbonds;
	frac1551 = (double)sum1551/sumbonds;
	frac1541 = (double)sum1541/sumbonds;
	frac1441 = (double)sum1441/sumbonds;
	frac1431 = (double)sum1431/sumbonds;
	frac142x = (double)sum142x/sumbonds;

	fp = fopen("pairCorrelationFunctions/CNA_information", "a+");
	fprintf(fp, "\nbond23\t");
	fprintf(fp, "%lf\t", frac1661);
	fprintf(fp, "%lf\t", frac1551);
	fprintf(fp, "%lf\t", frac1541);
	fprintf(fp, "%lf\t", frac1441);
	fprintf(fp, "%lf\t", frac1431);
	fprintf(fp, "%lf\t", frac142x);
	fclose(fp);
	//----------------------------------------------------------------

	return 0;
}
