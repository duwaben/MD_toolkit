#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>
#include <string>
#include <sys/stat.h>

using namespace std;

#define MAXBIN 1000 //the half of box length(multiply length by 2 in three dimensions)
#define delta_r 0.01
#define gap 1
#define Pi 3.141593
#define NOMIT 0 
#define NSTEP 30000
#define NATOM 200 //NATOM should not be set in my vision.

FILE *fp;
double LATTICE;
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];
double u1[NATOM], v1[NATOM], w1[NATOM];
double gr[MAXBIN];
double sumgr[MAXBIN];

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

void rearrangeCoordination(int bmax, int amin)
{
	double utemp[amin-bmax][NSTEP], vtemp[amin-bmax][NSTEP], wtemp[amin-bmax][NSTEP];
	for(int i = 0; i < NSTEP; ++ i)
	{
		for(int j = 0; j < amin-bmax; ++ j)
		{
			utemp[j][i] = u[j+bmax][i];
			vtemp[j][i] = v[j+bmax][i];
			wtemp[j][i] = w[j+bmax][i];
		}
	}
	
	for(int i = 0; i < NSTEP; ++ i)
	{
		for(int j = bmax; j < bmax+NATOM-amin; ++ j)
		{
			u[j][i] = u[j+amin-bmax][i];
			v[j][i] = v[j+amin-bmax][i];
			w[j][i] = w[j+amin-bmax][i];
		}
	}

	for(int i = 0; i < NSTEP; ++ i)
	{
		for(int j = NATOM-amin+bmax; j < NATOM; ++ j)
		{
			u[j][i] = utemp[j-NATOM+amin-bmax][i];
			v[j][i] = vtemp[j-NATOM+amin-bmax][i];
			w[j][i] = wtemp[j-NATOM+amin-bmax][i];
		}
	}

}

int main(int argc, char const *argv[])
{
	//namelist
	vector<string> namelist;
	readName(namelist);
	int step;
	//
	NumPerType thisNumPerType;
	thisNumPerType = readCoordinate("XDATCAR");
	mkdir("pairCorrelationFunctions", 00777);	

	int site = -1;
	string nameRe = "Re";
	for(int i = 0; i < namelist.size(); ++ i)
	{
		if(namelist[i] == nameRe)
		{
			site = i;
		}
	}

	if(site == -1)
	{
		cout << "No Re here." << endl;
		return 0;
	}

	if(site != namelist.size()-1 && site != 0)
	{
		rearrangeCoordination(thisNumPerType.getmin(site), thisNumPerType.getmax(site));
	}

	//X-X
	step = 0;
	for(int i = 0; i < MAXBIN; i ++)
		sumgr[i] = 0.0;

	for(int i = 0; i < NSTEP; i = i+gap){
		for(int j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		step ++;
		PDF_XX_8x(0, NATOM-thisNumPerType.getmax(site)+thisNumPerType.getmin(site));
		for(int j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j];
	}
	//step average
	for(int j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step;

	fp = fopen("pairCorrelationFunctions/PDF_X-X", "w");

	for(int k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * ((double)k+0.5), sumgr[k]);
	}
	fclose(fp);

	//Re-Re
	step = 0;
	for(int i = 0; i < MAXBIN; i ++)
		sumgr[i] = 0.0;

	for(int i = 0; i < NSTEP; i = i+gap){
		for(int j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		step ++;
		PDF_XX_8x(NATOM-thisNumPerType.getmax(site)+thisNumPerType.getmin(site), NATOM);
		for(int j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j];
	}
	//step average
	for(int j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step;

	fp = fopen("pairCorrelationFunctions/PDF_Re-Re", "w");

	for(int k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * ((double)k+0.5), sumgr[k]);
	}
	fclose(fp);

	//Re-X
	step = 0;
	for(int i = 0; i < MAXBIN; i ++)
		sumgr[i] = 0.0;

	for(int i = 0; i < NSTEP; i = i+gap){
		for(int j = 0; j < NATOM; j ++){
			u1[j] = u[j][i];
			v1[j] = v[j][i];
			w1[j] = w[j][i];
		}
		step ++;
		PDF_XY_8x(0, NATOM-thisNumPerType.getmax(site)+thisNumPerType.getmin(site), NATOM-thisNumPerType.getmax(site)+thisNumPerType.getmin(site), NATOM);
		for(int j = 0; j < MAXBIN; j ++)
			sumgr[j] += gr[j];
	}
	//step average
	for(int j = 0; j < MAXBIN; j ++)
			sumgr[j] /= step;

	fp = fopen("pairCorrelationFunctions/PDF_Re-X", "w");

	for(int k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * ((double)k+0.5), sumgr[k]);
	}
	fclose(fp);

	return 0;
}
