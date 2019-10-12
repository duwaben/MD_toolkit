#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NATOM 200
#define Lmin -11
#define Lmax 11

using namespace std;

int main(int argc, char const *argv[])
{
	FILE *fp;
	int NATOM_1, NATOM_2, NATOM_3;
	NATOM_1 = atoi(argv[1]);
	NATOM_2 = atoi(argv[2]);
	NATOM_3 = atoi(argv[3]);

	if(NATOM_1+NATOM_2+NATOM_3 != NATOM){
		cout << "Input error." << endl;
		exit(0);
	}

	struct DOS{
		double energy;
		int m;
		double s;
		double p;
		double d;
	};

	char name[7];
	char line[256];
	int n, i, j, k;
	struct DOS tmp, DOS1[(Lmax-Lmin)*50], DOS2[(Lmax-Lmin)*50], DOS3[(Lmax-Lmin)*50];
	double ENERGY[(Lmax-Lmin)*50], COUNT[(Lmax-Lmin)*50];
	int m[(Lmax-Lmin)*50];
	double energy, count;
	double Fermi;

	for(i = 0; i < (Lmax-Lmin)*50; i ++){
		ENERGY[i] = DOS1[i].energy = DOS2[i].energy = DOS3[i].energy = Lmin+i*0.02;
		DOS1[i].m = DOS2[i].m = DOS3[i].m = 0;
		DOS1[i].s = DOS1[i].p = DOS1[i].d = 0.0;
		DOS2[i].s = DOS2[i].p = DOS2[i].d = 0.0;
		DOS3[i].s = DOS3[i].p = DOS3[i].d = 0.0;
		COUNT[i] = 0.0;
		m[i] = 0;
	}


	for(n = 1; n <= 10; n ++){
		sprintf(name, "DOSCAR%d", n);
		fp = fopen(name, "r");

		for(i = 0; i < 5; ++ i)
			fgets(line, 256, fp);
		fscanf(fp, "%*lf%*lf%*d%lf", &Fermi);
		fgets(line, 256, fp);

		for(i = 0; i < 301; i ++){
			fscanf(fp, "%lf%lf%*lf", &energy, &count);
			fgets(line, 256, fp);
			energy -= Fermi;
			j = floor((energy-Lmin)/0.02);
			COUNT[j] += count;
			m[j] ++;
		}
		
		for(i = 0; i < 302; i ++)
			fgets(line, 256, fp);

		for(i = 0; i < NATOM_1; i ++){
			fgets(line, 512, fp);
			for(j = 0; j < 301; j ++){
				fscanf(fp, "%lf%lf%lf%lf", &tmp.energy, &tmp.s, &tmp.p, &tmp.d);
				fgets(line, 256, fp);
				tmp.energy -= Fermi;
                k = floor((tmp.energy-Lmin)/0.02);
				DOS1[k].m ++; DOS1[k].s += tmp.s; DOS1[k].p += tmp.p; DOS1[k].d += tmp.d; 
			}
		}

		for(i = 0; i < NATOM_2; i ++){
			fgets(line, 512, fp);
			for(j = 0; j < 301; j ++){
				fscanf(fp, "%lf%lf%lf%lf", &tmp.energy, &tmp.s, &tmp.p, &tmp.d);
				fgets(line, 256, fp);
				tmp.energy -= Fermi;
                k = floor((tmp.energy-Lmin)/0.02);
				DOS2[k].m ++; DOS2[k].s += tmp.s; DOS2[k].p += tmp.p; DOS2[k].d += tmp.d; 
			}
		}

		for(i = 0; i < NATOM_3; i ++){
			fgets(line, 512, fp);
			for(j = 0; j < 301; j ++){
				fscanf(fp, "%lf%lf%lf%lf", &tmp.energy, &tmp.s, &tmp.p, &tmp.d);
				fgets(line, 256, fp);
				tmp.energy -= Fermi;
                k = floor((tmp.energy-Lmin)/0.02);
				DOS3[k].m ++; DOS3[k].s += tmp.s; DOS3[k].p += tmp.p; DOS3[k].d += tmp.d; 
			}
		}

		fclose(fp);
	}

	for(i = 0; i < (Lmax-Lmin)*50; i ++){
		if(DOS1[i].m != 0){
			DOS1[i].s /= DOS1[i].m; DOS1[i].p /= DOS1[i].m; DOS1[i].d /= DOS1[i].m;
		}
		if(DOS2[i].m != 0){
			DOS2[i].s /= DOS2[i].m; DOS2[i].p /= DOS2[i].m; DOS2[i].d /= DOS2[i].m;
		}
		if(DOS3[i].m != 0){
			DOS3[i].s /= DOS3[i].m; DOS3[i].p /= DOS3[i].m; DOS3[i].d /= DOS3[i].m;
		}
		if(m[i] != 0)
			COUNT[i] /= m[i];
	}

//***************************************************************************************
//***************************************************************************************
	fp = fopen("DOS_tot", "w");
	for(i = 0; i < (Lmax-Lmin)*50; i ++){
		if(m[i] != 0)
		fprintf(fp, "%lf\t%lf\n", ENERGY[i], COUNT[i]);	
	}
	fclose(fp);
//---------------------------------------------------------------------------------------	
	fp = fopen("DOS_atom1", "w");
	for(i = 0; i < (Lmax-Lmin)*50; i ++){
		if(DOS1[i].m != 0)
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", DOS1[i].energy, DOS1[i].s, DOS1[i].p, DOS1[i].d);
	}
	fclose(fp);
//---------------------------------------------------------------------------------------	
	fp = fopen("DOS_atom2", "w");
	for(i = 0; i < (Lmax-Lmin)*50; i ++){
		if(DOS2[i].m != 0)
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", DOS2[i].energy, DOS2[i].s, DOS2[i].p, DOS2[i].d);
	}
	fclose(fp);
//---------------------------------------------------------------------------------------	
	fp = fopen("DOS_atom3", "w");
	for(i = 0; i < (Lmax-Lmin)*50; i ++){
		if(DOS3[i].m != 0)
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", DOS3[i].energy, DOS3[i].s, DOS3[i].p, DOS3[i].d);
	}
	fclose(fp);

	return 0;
}