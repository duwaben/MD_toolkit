/*
	Name: MSD(mean square displacement) and PAF(velocity autocorrelation functions) Processing Program
	Date: 28/03/18 
	Input: INCAR, XDATCAR, POTCAR;
	Output: MSD_1, MSD_2, MSD_3, DMSD_information, VAF_1, VAF_2, VAF_3, DVAF_information;
	Description: 
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>

#define TCOR 2000 //步数跨度（T0~T0+TCOR）
#define NOMIT 4000 //忽略步数
#define NSTEP 8000 //计算步数
#define NATOM 200
#define Pi  3.141593 

FILE*fp;
int TEBEG, POTIM;
int NATOM_1, NATOM_2, NATOM_3;
double LATTICE;
double NMASS_1, NMASS_2, NMASS_3;
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP]; 
double ru[NATOM][NSTEP], rv[NATOM][NSTEP], rw[NATOM][NSTEP];
double u1[NSTEP], v1[NSTEP], w1[NSTEP];
double MSD_1[TCOR], MSD_2[TCOR], MSD_3[TCOR], msddelta[TCOR];
double vu[NSTEP-1], vv[NSTEP-1], vw[NSTEP-1], vaf_a[TCOR], vaf_b[TCOR];

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

void velocityAutocorrelationFunction()
{
	int T, T0;
	double rafatemp, rafbtemp;
	int TMAX;
	for(T = 0; T < TCOR; T ++){
		rafatemp = rafbtemp = 0.0;
		TMAX = NSTEP - 1 - T;
		for(T0 = 0; T0 < TMAX; T0 ++){
			rafatemp += vu[T0] * vu[T0+T] + vv[T0] * vv[T0+T] + vw[T0] * vw[T0+T];
			rafbtemp += vu[T0] * vu[T0] + vv[T0] * vv[T0] + vw[T0] * vw[T0];
		}
		vaf_a[T] = rafatemp/(double)TMAX; 
		vaf_b[T] = rafbtemp/(double)TMAX;
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

void readINCAR(const char *name)
{
	char *word1 = "POTIM";
	char *word2 = "=";
	char *word3 = "TEBEG";

	fp = fopen(name, "r");
	findWord(word1);
	findWord(word2);
	fscanf(fp, "%d", &POTIM);
	fclose(fp);

	fp = fopen(name, "r");
	findWord(word3);
	findWord(word2);
	fscanf(fp, "%d", &TEBEG);
	fclose(fp);
}

void readPOTCAR(const char *name)
{
	fp = fopen(name, "r");
	const char *word1 = "POMASS";
	const char *word2 = "=";
	findWord(word1);findWord(word2);
	fscanf(fp, "%lf", &NMASS_1);
	findWord(word1);findWord(word2);
	fscanf(fp, "%lf", &NMASS_2);
	findWord(word1);findWord(word2);
	fscanf(fp, "%lf", &NMASS_3);
	fclose(fp);
}

int main()
{  
	int i, j, T;
			
	readCoordinate("XDATCAR");
	readINCAR("INCAR");
	readPOTCAR("POTCAR");
	mkdir("dynamicalProperty", 00777);	

	//坐标修正：用于计算MSD和VAF
	for(i = 0; i < NSTEP; i ++){
		for( j = 0; j < NATOM; j ++){
			ru[j][i] = u[j][i];
			rw[j][i] = w[j][i];
			rv[j][i] = v[j][i];
		}
	}

	for(i = 0; i < NSTEP - 1; i ++){
		for( j = 0; j < NATOM; j ++){
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

	//1st type
	for(j = 0; j < NATOM_1; j ++){
		for(i = 0; i < NSTEP; i ++){
			u1[i] = ru[j][i];	
			v1[i] = rv[j][i];	
			w1[i] = rw[j][i];
		}
		meanSquareDisplacement(); 
		for(T = 0; T < TCOR; T ++){
			MSD_1[T] += msddelta[T];
		}
	}

	for(T = 0; T < TCOR; T ++){
		MSD_1[T] = MSD_1[T]/NATOM_1;
	}
	
	//2ed type
	for(j = NATOM_1; j < NATOM_1+NATOM_2; j ++){
		for(i = 0; i < NSTEP; i ++){
			u1[i] = ru[j][i];	
			v1[i] = rv[j][i];	
			w1[i] = rw[j][i];
		}
		meanSquareDisplacement(); 
		for(T = 0; T < TCOR; T ++){
			MSD_2[T] += msddelta[T];
		}
	}
	for(T = 0; T < TCOR; T ++){
		MSD_2[T] = MSD_2[T]/NATOM_2;
	}
 
	//3rd type
	for(j = NATOM_1+NATOM_2; j < NATOM; j ++){
		for(i = 0; i < NSTEP; i ++){
			u1[i] = ru[j][i];	
			v1[i] = rv[j][i];	
			w1[i] = rw[j][i];
		}
		meanSquareDisplacement(); 
		for(T = 0; T < TCOR; T ++){
			MSD_3[T] += msddelta[T];
		}
	}
	for(T = 0; T < TCOR; T ++){
		MSD_3[T] = MSD_3[T]/NATOM_3;
	}

	fp = fopen("dynamicalProperty/MSD_1", "w");
	for(T = 0; T < TCOR; T++){
		fprintf(fp, "%d\t%lf\n", T, MSD_1[T]);
	}
	fclose(fp);
	
	fp = fopen("dynamicalProperty/MSD_2", "w");
	for(T = 0; T < TCOR; T ++){
		fprintf(fp, "%d\t%lf\n", T, MSD_2[T]);
	}
	fclose(fp);
	
	fp = fopen("dynamicalProperty/MSD_3", "w");
	for(T = 0; T < TCOR; T++){
		fprintf(fp, "%d\t%lf\n", T, MSD_3[T]);
	}
	fclose(fp);


	double DMSD_1, DMSD_2, DMSD_3;
	double area_1, area_2, area_3;
	area_1 = area_2 = area_3 = 0.0;
	for(i = 0; i < TCOR; i ++){
		area_1 += 0.001 * POTIM * MSD_1[i];
		area_2 += 0.001 * POTIM * MSD_2[i];
		area_3 += 0.001 * POTIM * MSD_3[i];
	}
	double bottomEdge = TCOR * 0.001 * POTIM;
	DMSD_1 = area_1 * 2 /(bottomEdge * bottomEdge) /6;
	DMSD_2 = area_2 * 2 /(bottomEdge * bottomEdge) /6;
	DMSD_3 = area_3 * 2 /(bottomEdge * bottomEdge) /6;


	fp = fopen("dynamicalProperty/DMSD_information", "w");
	fprintf(fp, "atom1: %lf\n", DMSD_1);
	fprintf(fp, "atom2: %lf\n", DMSD_2);
	fprintf(fp, "atom3: %lf\n", DMSD_3);
	fclose(fp);
	


	//velocity autocorrelation function
	double VAFA[TCOR], VAFB[TCOR], quotient_1[TCOR], quotient_2[TCOR], quotient_3[TCOR];
	double integration_1, integration_2, integration_3;

	for(i = 0; i < TCOR; ++ i){
		VAFA[i] = VAFB[i] = 0;
	}
	for(j = 0; j < NATOM_1; j ++){
		for(i = 0; i < NSTEP-1; i ++){
			vu[i] = ru[j][i+1] - ru[j][i];
			vu[i] /= 0.001 * POTIM;
			vv[i] = rv[j][i+1] - rv[j][i];
			vv[i] /= 0.001 * POTIM;
			vw[i] = rw[j][i+1] - rw[j][i];
			vw[i] /= 0.001 * POTIM;
		} 
		velocityAutocorrelationFunction();
		for(T = 0; T < TCOR; T ++){
			VAFA[T] += vaf_a[T];
			VAFB[T] += vaf_b[T];
		}	
	}
	for(T = 0; T < TCOR; T ++){
			// VAFA[T] /= NATOM_1;
			// VAFB[T] /= NATOM_1;
			quotient_1[T] = VAFA[T]/VAFB[T];
	}


	for(i = 0; i < TCOR; ++ i){
		VAFA[i] = VAFB[i] = 0;
	}
	for(j = NATOM_1; j < NATOM_1+NATOM_2; j ++){
		for(i = 0; i < NSTEP-1; i ++){
			vu[i] = ru[j][i+1] - ru[j][i];
			vu[i] /= 0.001 * POTIM;
			vv[i] = rv[j][i+1] - rv[j][i];
			vv[i] /= 0.001 * POTIM;
			vw[i] = rw[j][i+1] - rw[j][i];
			vw[i] /= 0.001 * POTIM;
		} 
		velocityAutocorrelationFunction();
		for(T = 0; T < TCOR; T ++){
			VAFA[T] += vaf_a[T];
			VAFB[T] += vaf_b[T];
		}	
	}
	for(T = 0; T < TCOR; T ++){
			// VAFA[T] /= NATOM_2;
			// VAFB[T] /= NATOM_2;
			quotient_2[T] = VAFA[T]/VAFB[T];
	}
	

	for(i = 0; i < TCOR; ++ i){
		VAFA[i] = VAFB[i] = 0;
	}
	for(j = NATOM_1+NATOM_2; j < NATOM; j ++){
		for(i = 0; i < NSTEP-1; i ++){
			vu[i] = ru[j][i+1] - ru[j][i];
			vu[i] /= 0.001 * POTIM;
			vv[i] = rv[j][i+1] - rv[j][i];
			vv[i] /= 0.001 * POTIM;
			vw[i] = rw[j][i+1] - rw[j][i];
			vw[i] /= 0.001 * POTIM;
		} 
		velocityAutocorrelationFunction();
		for(T = 0; T < TCOR; T ++){
			VAFA[T] += vaf_a[T];
			VAFB[T] += vaf_b[T];
		}	
	}
	for(T = 0; T < TCOR; T ++){
			// VAFA[T] /= NATOM_3;
			// VAFB[T] /= NATOM_3;
			quotient_3[T] = VAFA[T]/VAFB[T];
	}


	fp = fopen("dynamicalProperty/VAF_1", "w");
	for(T = 0; T < TCOR; T ++){
		fprintf(fp, "%d\t%lf\n", T, quotient_1[T]);
	}
	fclose(fp);
	
	fp = fopen("dynamicalProperty/VAF_2", "w");
	for(T = 0; T < TCOR; T ++){
		fprintf(fp, "%d\t%lf\n", T, quotient_2[T]);
	}
	fclose(fp);

	fp = fopen("dynamicalProperty/VAF_3", "w");
	for(T = 0; T < TCOR; T ++){
		fprintf(fp, "%d\t%lf\n", T, quotient_3[T]);
	}
	fclose(fp);


	integration_1 = integration_2 = integration_3 = 0.0;
	for(T = 0; T < TCOR; T ++){
		integration_1 += quotient_1[T]; integration_2 += quotient_2[T]; integration_3 += quotient_3[T];
	}
	integration_1 *= 0.8282 * (0.001 * POTIM) * TEBEG/NMASS_1;
	integration_2 *= 0.8282 * (0.001 * POTIM) * TEBEG/NMASS_2;
	integration_3 *= 0.8282 * (0.001 * POTIM) * TEBEG/NMASS_3;

	fp=fopen("dynamicalProperty/DVAF_information", "w");
	fprintf(fp, "atom1: %lf\n", integration_1);
	fprintf(fp, "atom2: %lf\n", integration_2);
	fprintf(fp, "atom3: %lf", integration_3);
	fclose(fp);

	return 0;	
}



