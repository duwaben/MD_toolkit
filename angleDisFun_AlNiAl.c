#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NATOM 200       //total atom number
#define NOMIT 4000 		//the number of omitted steps
#define NSTEP 8000		//the number of counting steps
#define gap 1			//the interval of counting
#define delt_a 0.5
#define MAXBIN 360
#define Pi 3.141592653

FILE *fp;
double LATTICE;
int NATOM_1, NATOM_2, NATOM_3; //three types
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];	
double u1[NATOM], v1[NATOM], w1[NATOM];

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

double bondAngle(int atom1, int atom2, int atom3, double cutoff1, double cutoff2)
{
	double r12 = r(u1[atom1], v1[atom1], w1[atom1], u1[atom2], v1[atom2], w1[atom2]);
	double r13 = r(u1[atom1], v1[atom1], w1[atom1], u1[atom3], v1[atom3], w1[atom3]);
	if( r12 > cutoff1 || r13 > cutoff2 || r12 < 0.001 || r13 < 0.001 )
		return -1.0;
	double r23 = r(u1[atom2], v1[atom2], w1[atom2], u1[atom3], v1[atom3], w1[atom3]);
	double anglecos = ( r12 * r12 + r13 * r13 - r23 * r23 ) /( 2 * r12 * r13 );
	return acos(anglecos)*180/Pi;
}

int main(int argc, char *argv[])
{	

	int i, j, k;
	readCoordinate("XDATCAR");

	double angle;
	int sumg_angle;
	int g_angle[MAXBIN];
	double g_angleDistribution[MAXBIN];
  	int BIN;
	int istep;

	double CUTOFF;
	CUTOFF = atof(argv[1]);
//---------------------------------------------------------------------------
	//Ni(Al,Al)
	for( BIN = 0; BIN < MAXBIN; BIN ++ ){
		g_angle[BIN] = 0;
	}
	sumg_angle = 0;

	for( istep = 0; istep < NSTEP; istep += gap ){
		
		for( i = 0; i < NATOM; i ++ ){
			u1[i] = u[i][istep];
			v1[i] = v[i][istep];
			w1[i] = w[i][istep];
		}

		for( i = 0; i < NATOM_1; i ++ )
			for( j = NATOM_1; j < NATOM_1+NATOM_2; j ++ )
				for( k = j+1; k < NATOM_1+NATOM_2; k ++ ){
					
					angle = bondAngle( i, j, k, CUTOFF, CUTOFF );
					if( angle < 0 )
						continue;
					
					BIN = floor(angle /delt_a);
					g_angle[BIN] ++;				
				}
	}

	for( BIN = 0; BIN < MAXBIN; BIN ++ ){
		sumg_angle += g_angle[BIN];
	}
	for( BIN = 0; BIN < MAXBIN; BIN ++ ){
		g_angleDistribution[BIN] = (double)g_angle[BIN] /(double)sumg_angle /delt_a;
	}

	fp = fopen( "angleDisFun_Ni_AlAl", "w" );
	for( BIN = 0; BIN < MAXBIN; BIN ++){
		fprintf(fp, "%lf\t%lf\n", (BIN+0.5)*delt_a, g_angleDistribution[BIN]);
	}
	fclose(fp);

//------------------------------------------------------------------------------	

	return 0;
}
