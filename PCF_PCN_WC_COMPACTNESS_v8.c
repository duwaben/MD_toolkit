/*
	Name: Pair Correlation Function calculation program 8.1
	Time: 20/03/2018
	Input file: XDATCAR
	Output file: PDF11_8x, PDF22_8x, PDF33_8x, PDF12_8x, PDF13_8x, PDF23_8x, 
				CUTOFFPEAK_information, pCN_information, WC_information, COMPACTNESS_information.
	Improvement: read x coordinate of first peak and first vale automatically.
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
#define NOMIT 0 
#define NSTEP 30000
#define NATOM 200

FILE *fp;
double LATTICE;
int NATOM_1, NATOM_2, NATOM_3;
double CUTOFF11, CUTOFF12, CUTOFF13, CUTOFF22, CUTOFF23, CUTOFF33;
double PEAK11, PEAK12, PEAK13, PEAK22, PEAK23, PEAK33;
double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];
double u1[NATOM], v1[NATOM], w1[NATOM];

double gr[MAXBIN];
double sumgr[MAXBIN];
double fft_filter[MAXBIN];

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

int gaussSolve(a, b)
double a[], b[];
{
	int *js, l, k, i, j, is, p, q;
    double d, t;
    js = malloc(4 * sizeof(int));
    l = 1;
    for(k = 0; k <= 2; k ++){ 
    	d = 0.0;
        for(i = k; i <= 3; i ++)
          	for(j = k; j <= 3; j ++){ 
          		t = fabs(a[i*4+j]);
              	if(t > d){ 
              		d = t; 
              		js[k] = j; 
              		is = i;
              	}
            }
        if(d + 1.0 == 1.0) 
        	l = 0;
        else{ 
        	if(js[k] != k)
              	for(i = 0; i <= 3; i ++){ 
              		p = i*4 + k; q = i*4 + js[k];
                  	t = a[p]; a[p] = a[q]; a[q] = t;
                }
            if(is != k){ 
            	for(j = k; j <= 3; j ++){ 
            		p = k*4 + j; q = is*4 + j;
                    t = a[p]; a[p] = a[q]; a[q] = t;
                }
                t = b[k]; b[k] = b[is]; b[is] = t;
            }
        }
        if(l == 0){ 
        	free(js); 
        	printf("fail\4");
            return(0);
        }
        d = a[k*4+k];
        for(j = k+1; j <= 3; j ++){ 
        	p = k*4 + j; a[p] = a[p]/d;
        }
        b[k] = b[k]/d;
        for(i = k+1; i <= 3; i ++){ 
        	for(j = k+1; j <= 3; j ++){ 
        		p = i*4 + j;
                a[p] = a[p] - a[i*4+k]*a[k*4+j];
            }
            b[i] = b[i] - a[i*4+k]*b[k];
        }
    }
    d = a[15];
    if(fabs(d) + 1.0 == 1.0){ 
    	free(js); printf("fail\4");
        return(0);
    }
    b[3] = b[3] /d;
    for(i = 2; i >= 0; i --){ 
    	t = 0.0;
        for(j = i+1; j <= 3; j ++)
          	t = t + a[i*4+j]*b[j];
        b[i] = b[i]-t;
    }
    js[3] = 3;
    for(k = 3; k >= 0; k --)
      	if(js[k] != k){ 
      		t = b[k]; b[k] = b[js[k]]; b[js[k]] = t;
      	}
    free(js);
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

double peakParameter(int k)
{
	double ans;
	if(sumgr[k+25] > sumgr[k])
		ans = sumgr[k+25] - sumgr[k];
	else
		ans = sumgr[k] - sumgr[k+25];
	return ans;
}

double fitCubicCurve(int x2, int peak)
{
	// 2+2
	int x1 = x2 - 10;
	int x3 = x2 + 25;
	int x4 = x3 + 10;
	double coordinate_peak;
	double A[4][4], array[4];
	A[0][0] = x1 * x1 * x1; A[0][1] = x1 * x1; A[0][2] = x1; A[0][3] = 1;
	A[1][0] = x2 * x2 * x2; A[1][1] = x2 * x2; A[1][2] = x2; A[1][3] = 1;
	A[2][0] = x3 * x3 * x3; A[2][1] = x3 * x3; A[2][2] = x3; A[2][3] = 1;
	A[3][0] = x4 * x4 * x4; A[3][1] = x4 * x4; A[3][2] = x4; A[3][3] = 1;
	array[0] = sumgr[x1]; 
	array[1] = sumgr[x2]; 
	array[2] = sumgr[x3]; 
	array[3] = sumgr[x4]; 

	gaussSolve(A, array);
	if(peak)
		coordinate_peak = (- array[1] - sqrt(array[1] * array[1] - 3 * array[0] * array[2])) /(3 * array[0]);
	else
		coordinate_peak = (- array[1] + sqrt(array[1] * array[1] - 3 * array[0] * array[2])) /(3 * array[0]);

	return (coordinate_peak + 0.5) * delta_r;		
}

void fftFiltrate()
{
	int i, j;
	double ReDFT_out[MAXBIN/2+1], ImDFT_out[MAXBIN/2+1];
	double sum1, sum2;
	for(i = 0; i < MAXBIN/2+1; ++ i){
		sum1 = 0.0;
		sum2 = 0.0;
		for(j = 0; j < MAXBIN; j ++){
			sum1 += sumgr[j] * cos(2 * Pi * i * j /MAXBIN); 
			sum2 += sumgr[j] * sin(2 * Pi * i * j /MAXBIN); 
		}
		ReDFT_out[i] = 2.0/MAXBIN * sum1;
		ImDFT_out[i] = -2.0/MAXBIN * sum2;
	}
	ReDFT_out[0] /= 2.0;
	ImDFT_out[0] /= -2.0;	


	for(i = 40; i < MAXBIN/2+1; ++ i){
		ReDFT_out[i] = 0.0;
		ImDFT_out[i] = 0.0;
	}

	double sum;
	for(i = 0; i < MAXBIN; i ++){
		sum = 0.0;
		for(j = 0; j < MAXBIN/2+1; j ++){
			sum += ReDFT_out[j] * cos(2 * Pi * j * i /MAXBIN) - ImDFT_out[j] * sin(2 * Pi * j * i /MAXBIN);
		}
		fft_filter[i] = sum;
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

int main()
{
	int i, j, k;
	double pCN_11, pCN_12, pCN_13, pCN_21, pCN_22, pCN_23, pCN_31, pCN_32, pCN_33;
	readCoordinate("XDATCAR");
	mkdir("pairCorrelationFunctions", 00777);	
 //	mkdir("pairCorrelationFunctions", 0777);	in windows OS





	//*********************************************************
	//*********************************************************
	
	//pair_11
	int step = 0; 
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

	if((fp = fopen("pairCorrelationFunctions/PDF11_8x","w")) == NULL){
		printf("Cannot create file PDF11."); exit(0);
	}

	for(k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * ((double)k+0.5), sumgr[k]);
	}
	fclose(fp);

	//find peak and vale
	//-------------------------------------------------------------------------
	int peakflag = 1;
	int	xpeak = 0;
	while(peakflag && xpeak < MAXBIN){
		peakflag = 3;
		for(j = xpeak; j < xpeak+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] < sumgr[j])
				peakflag --;
		}
		xpeak ++;
	}
	xpeak --;
	peakflag = 1;
	int xvale = xpeak + 1;
	while(peakflag && xvale < MAXBIN){
		peakflag = 3;
		for(j = xvale; j < xvale+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] > sumgr[j])
				peakflag --;
		}
		xvale ++;
	}
	xvale --;

	fftFiltrate();
	int xvale_fft;
	double cutoff_fft;
	int xvale_fft_l;
	int xvale_fft_r;
	xvale_fft = xvale;
	for(i = xvale; i < xvale+25; i ++){
		if(fft_filter[i] < fft_filter[xvale_fft])
			xvale_fft = i;
	}
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i --;
	xvale_fft_l = i;
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i ++;
	xvale_fft_r = i;
	cutoff_fft = delta_r * ((xvale_fft_l + xvale_fft_r)/2.0 + 0.5);

	double temp;
	PEAK11 = fitCubicCurve(xpeak, 1);
	temp = fitCubicCurve(xvale, 0);
	if(temp - cutoff_fft > 0.07 || cutoff_fft - temp > 0.07)
		CUTOFF11 = cutoff_fft;
	else
		CUTOFF11 = temp;
	//--------------------------------------------------------------------------

	//pCN
	double sumpCN = 0;
	int CUTOFFBIN = floor(CUTOFF11/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_11 = sumpCN * 4 * Pi * NATOM_1 /(LATTICE * LATTICE * LATTICE);
	//*********************************************************
	//*********************************************************





	//pair_22
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

	if((fp = fopen("pairCorrelationFunctions/PDF22_8x","w")) == NULL){
		printf("Cannot create file PDF22."); exit(0);
	}

	for(k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * ((double)k+0.5), sumgr[k]);
	}
	fclose(fp);

	//find peak and vale
	//-------------------------------------------------------------------------
	peakflag = 1;
	xpeak = 0;
	while(peakflag && xpeak < MAXBIN){
		peakflag = 3;
		for(j = xpeak; j < xpeak+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] < sumgr[j])
				peakflag --;
		}
		xpeak ++;
	}
	xpeak --;
	peakflag = 1;
	xvale = xpeak + 1;
	while(peakflag && xvale < MAXBIN){
		peakflag = 3;
		for(j = xvale; j < xvale+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] > sumgr[j])
				peakflag --;
		}
		xvale ++;
	}
	xvale --;

	fftFiltrate();
	xvale_fft = xvale;
	for(i = xvale; i < xvale+25; i ++){
		if(fft_filter[i] < fft_filter[xvale_fft])
			xvale_fft = i;
	}
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i --;
	xvale_fft_l = i;
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i ++;
	xvale_fft_r = i;
	cutoff_fft = delta_r * ((xvale_fft_l + xvale_fft_r)/2.0 + 0.5);

	PEAK22 = fitCubicCurve(xpeak, 1);
	temp = fitCubicCurve(xvale, 0);
	if(temp - cutoff_fft > 0.07 || cutoff_fft - temp > 0.07)
		CUTOFF22 = cutoff_fft;
	else
		CUTOFF22 = temp;
	//--------------------------------------------------------------------------

	//pCN
	sumpCN = 0;
	CUTOFFBIN = floor(CUTOFF22/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_22 = sumpCN * 4 * Pi * NATOM_2 /(LATTICE * LATTICE * LATTICE);
	//*********************************************************
	//*********************************************************





	//pair_33
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

	if((fp = fopen("pairCorrelationFunctions/PDF33_8x","w")) == NULL){
		printf("Cannot create file PDF33."); exit(0);
	}

	for(k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * ((double)k+0.5), sumgr[k]);
	}
	fclose(fp);

	//find peak and vale
	//-------------------------------------------------------------------------
	peakflag = 1;
	xpeak = 0;
	while(peakflag && xpeak < MAXBIN){
		peakflag = 3;
		for(j = xpeak; j < xpeak+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] < sumgr[j])
				peakflag --;
		}
		xpeak ++;
	}
	xpeak --;
	peakflag = 1;
	xvale = xpeak + 1;
	while(peakflag && xvale < MAXBIN){
		peakflag = 3;
		for(j = xvale; j < xvale+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] > sumgr[j])
				peakflag --;
		}
		xvale ++;
	}
	xvale --;

	fftFiltrate();
	xvale_fft = xvale;
	for(i = xvale; i < xvale+25; i ++){
		if(fft_filter[i] < fft_filter[xvale_fft])
			xvale_fft = i;
	}
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i --;
	xvale_fft_l = i;
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i ++;
	xvale_fft_r = i;
	cutoff_fft = delta_r * ((xvale_fft_l + xvale_fft_r)/2.0 + 0.5);

	PEAK33 = fitCubicCurve(xpeak, 1);
	temp = fitCubicCurve(xvale, 0);
	if(temp - cutoff_fft > 0.07 || cutoff_fft - temp > 0.07)
		CUTOFF33 = cutoff_fft;
	else
		CUTOFF33 = temp;
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------

	//pCN
	sumpCN = 0;
	CUTOFFBIN = floor(CUTOFF33/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_33 = sumpCN * 4 * Pi * NATOM_3 /(LATTICE * LATTICE * LATTICE);
	//*********************************************************
	//*********************************************************





	//PDF_XY
	//pair_12
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

	if((fp = fopen("pairCorrelationFunctions/PDF12_8x","w")) == NULL){
		printf("Cannot create file PDF12."); exit(0);
	}

	for(k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * ((double)k+0.5), sumgr[k]);
	}
	fclose(fp);

	//find peak and vale
	//-------------------------------------------------------------------------
	peakflag = 1;
	xpeak = 0;
	while(peakflag && xpeak < MAXBIN){
		peakflag = 3;
		for(j = xpeak; j < xpeak+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] < sumgr[j])
				peakflag --;
		}
		xpeak ++;
	}
	xpeak --;
	peakflag = 1;
	xvale = xpeak + 1;
	while(peakflag && xvale < MAXBIN){
		peakflag = 3;
		for(j = xvale; j < xvale+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] > sumgr[j])
				peakflag --;
		}
		xvale ++;
	}
	xvale --;

	fftFiltrate();
	xvale_fft = xvale;
	for(i = xvale; i < xvale+25; i ++){
		if(fft_filter[i] < fft_filter[xvale_fft])
			xvale_fft = i;
	}
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i --;
	xvale_fft_l = i;
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i ++;
	xvale_fft_r = i;
	cutoff_fft = delta_r * ((xvale_fft_l + xvale_fft_r)/2.0 + 0.5);

	PEAK12 = fitCubicCurve(xpeak, 1);
	temp = fitCubicCurve(xvale, 0);
	if(temp - cutoff_fft > 0.07 || cutoff_fft - temp > 0.07)
		CUTOFF12 = cutoff_fft;
	else
		CUTOFF12 = temp;
	//--------------------------------------------------------------------------

	//pCN
	sumpCN = 0;
	CUTOFFBIN = floor(CUTOFF12/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_12 = sumpCN * 4 * Pi * NATOM_2 /(LATTICE * LATTICE * LATTICE);
	pCN_21 = sumpCN * 4 * Pi * NATOM_1 /(LATTICE * LATTICE * LATTICE);

	//*********************************************************
	//*********************************************************





	//pair_13
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

	if((fp = fopen("pairCorrelationFunctions/PDF13_8x","w")) == NULL){
		printf("Cannot create file PDF13."); exit(0);
	}

	for(k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * ((double)k+0.5), sumgr[k]);
	}
	fclose(fp);

	//find peak and vale
	//-------------------------------------------------------------------------
	peakflag = 1;
	xpeak = 0;
	while(peakflag && xpeak < MAXBIN){
		peakflag = 3;
		for(j = xpeak; j < xpeak+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] < sumgr[j])
				peakflag --;
		}
		xpeak ++;
	}
	xpeak --;
	peakflag = 1;
	xvale = xpeak + 1;
	while(peakflag && xvale < MAXBIN){
		peakflag = 3;
		for(j = xvale; j < xvale+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] > sumgr[j])
				peakflag --;
		}
		xvale ++;
	}
	xvale --;

	fftFiltrate();
	xvale_fft = xvale;
	for(i = xvale; i < xvale+25; i ++){
		if(fft_filter[i] < fft_filter[xvale_fft])
			xvale_fft = i;
	}
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i --;
	xvale_fft_l = i;
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i ++;
	xvale_fft_r = i;
	cutoff_fft = delta_r * ((xvale_fft_l + xvale_fft_r)/2.0 + 0.5);

	PEAK13 = fitCubicCurve(xpeak, 1);
	temp = fitCubicCurve(xvale, 0);
	if(temp - cutoff_fft > 0.07 || cutoff_fft - temp > 0.07)
		CUTOFF13 = cutoff_fft;
	else
		CUTOFF13 = temp;

	//pCN
	sumpCN = 0;
	CUTOFFBIN = floor(CUTOFF13/delta_r);
	for(k = 0; k < CUTOFFBIN; k ++){
		sumpCN += (delta_r * (k+0.5)) * (delta_r * (k+0.5)) * sumgr[k] * delta_r;
	}
	pCN_13 = sumpCN * 4 * Pi * NATOM_3 /(LATTICE * LATTICE * LATTICE);
	pCN_31 = sumpCN * 4 * Pi * NATOM_1 /(LATTICE * LATTICE * LATTICE);

	//*********************************************************
	//*********************************************************





	//pair_23

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

	if((fp = fopen("pairCorrelationFunctions/PDF23_8x","w")) == NULL){
		printf("Cannot create file PDF23."); exit(0);
	}

	for(k = 0; k < MAXBIN; k ++){
		fprintf(fp, "%lf\t%lf\n", delta_r * ((double)k+0.5), sumgr[k]);
	}
	fclose(fp);

	//find peak and vale
	//-------------------------------------------------------------------------
	peakflag = 1;
	xpeak = 0;
	while(peakflag && xpeak < MAXBIN){
		peakflag = 3;
		for(j = xpeak; j < xpeak+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] < sumgr[j])
				peakflag --;
		}
		xpeak ++;
	}
	xpeak --;
	peakflag = 1;
	xvale = xpeak + 1;
	while(peakflag && xvale < MAXBIN){
		peakflag = 3;
		for(j = xvale; j < xvale+5; j ++){
			if(peakParameter(j) < peakParameter(j+1) && sumgr[j+25] > sumgr[j])
				peakflag --;
		}
		xvale ++;
	}
	xvale --;

	fftFiltrate();
	xvale_fft = xvale;
	for(i = xvale; i < xvale+25; i ++){
		if(fft_filter[i] < fft_filter[xvale_fft])
			xvale_fft = i;
	}
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i --;
	xvale_fft_l = i;
	i = xvale_fft;
	while(fft_filter[i] < fft_filter[xvale_fft] + 0.01)
		i ++;
	xvale_fft_r = i;
	cutoff_fft = delta_r * ((xvale_fft_l + xvale_fft_r)/2.0 + 0.5);

	PEAK23 = fitCubicCurve(xpeak, 1);
	temp = fitCubicCurve(xvale, 0);
	if(temp - cutoff_fft > 0.07 || cutoff_fft - temp > 0.07)
		CUTOFF23 = cutoff_fft;
	else
		CUTOFF23 = temp;
	//--------------------------------------------------------------------------

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

	//***********************************************************************
	//***********************************************************************
	//record the information of CUTOFF and coordinate of first peak.
	fp = fopen("pairCorrelationFunctions/CUTOFFPEAK_information", "w");

	fprintf(fp, "CUTOFF11: %lf\n", CUTOFF11);
	fprintf(fp, "CUTOFF12: %lf\n", CUTOFF12);
	fprintf(fp, "CUTOFF13: %lf\n", CUTOFF13);
	fprintf(fp, "CUTOFF22: %lf\n", CUTOFF22);
	fprintf(fp, "CUTOFF23: %lf\n", CUTOFF23);
	fprintf(fp, "CUTOFF33: %lf\n", CUTOFF33);

	fprintf(fp, "PEAK11: %lf\n", PEAK11);
	fprintf(fp, "PEAK12: %lf\n", PEAK12);
	fprintf(fp, "PEAK13: %lf\n", PEAK13);
	fprintf(fp, "PEAK22: %lf\n", PEAK22);
	fprintf(fp, "PEAK23: %lf\n", PEAK23);
	fprintf(fp, "PEAK33: %lf\n", PEAK33);

	fclose(fp);
	//***********************************************************************
	//***********************************************************************


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
