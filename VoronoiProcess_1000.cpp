/*
	Name: Voronoi Processing Program_1000
	Date: 17/04/18 
	Input: XDATCAR;
	Output: voroimport, RESULT, VORONOI_RESULT;
	Description: improve the criterion to remove small surfaces
*/

//it is noted that the last parameter in container constructor should be carefully set
#include "voro++.hh"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

using namespace voro;

#define NOMIT 9000
#define NSTEP 3000
#define NATOM 200
// #define AreaLimit 0.001

void findWord(char * word, FILE * fp)
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

void voronoiGenerate()
{
	FILE * fp;
	container * con;

	char line[1024];
	char neighbourstring[512];
	char areastring[512];
	int i, j, k, t;
	int numline;
	double pos[3];

	const double x_min = 0, x_max = 1;
	const double y_min = 0, y_max = 1;
	const double z_min = 0, z_max = 1;
	const int n_x = 3, n_y = 3, n_z = 3;

	struct VOROPOLY{
		int id; double pos[3]; int count[11]; int id_n[25]; double area[25]; double Area_total;
	}; 

	VOROPOLY * voro;
	voro = new VOROPOLY[NATOM];
	for(i = 0; i < NATOM; i ++){
		for(j = 0; j < 11; j ++)
			voro[i].count[j] = -1;
		for(j = 0; j < 25; j ++)
			voro[i].id_n[j] = -1;
		for(j = 0; j < 25; j ++)
			voro[i].area[j] = -1.0;
	}

	if(fopen("voronoiTessellation/RESULT", "r") != NULL)
		remove("voronoiTessellation/RESULT");
	
	con = new container(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, true, true, true, 8);	
	con->import("voronoiTessellation/voroimport");
	con->print_custom("ID =%i\t ||Position =%x %y %z\t||ID_neighbour =%n\t||Area_face =%f\t||Area_total =%F\t||", "voronoiTessellation/result1");
	delete con; 

	fp = fopen("voronoiTessellation/result1","r");
    for(i = 0; i < NATOM; i ++){
    	fgets(line, 1024, fp);
    	sscanf(line, "%*[^=]=%d%*[^=]=%lf %lf %lf%*[^=]=%[^|]%*[^=]=%[^|]%*[^=]=%lf", &voro[i].id, &voro[i].pos[0], &voro[i].pos[1], &voro[i].pos[2], neighbourstring, areastring, &voro[i].Area_total);
		
		k = 0;
		for(j = 0; neighbourstring[j]; j ++){
			if(k >= 25)printf("Warning: overflow\n");
			sscanf(neighbourstring+j, "%d%n", &voro[i].id_n[k++], &t); j += t;
		}

		k = 0;
		for(j = 0; areastring[j]; j ++){
			if(k >= 25)printf("Warning: overflow\n");
			sscanf(areastring+j, "%lf%n", &voro[i].area[k++], &t); j += t;
		}
	}
	fclose(fp);
	for(i = 0; i < NATOM; i ++){
		numline = 1;
		con = new container(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, true, true, true, 8);
		for(j = 0; j < 25; j ++){
			for(k = 0; k < NATOM; k ++){
				if(voro[i].id_n[j] == voro[k].id){
					pos[0] = voro[k].pos[0]; pos[1] = voro[k].pos[1]; pos[2] = voro[k].pos[2];
				}
			}			
			if(voro[i].area[j] > voro[i].Area_total * 0.01){
				con->put(voro[i].id_n[j], pos[0], pos[1], pos[2]); numline++;
			}
		}

		con->put(voro[i].id, voro[i].pos[0], voro[i].pos[1], voro[i].pos[2]);
		con->print_custom("ID =%i\t||Freq_table =%A\t||ID_neighbour =%n\t||volume =%v", "voronoiTessellation/result2");
		delete con; con = NULL;
        
		fp = fopen("voronoiTessellation/result2","r");
		for(j = 0; j < numline; j ++){
			fgets(line, 1024, fp);
			sscanf(line, "%*[^=]=%d", &k);
			if(k == voro[i].id)break;
		}
		fclose(fp);

		fp = fopen("voronoiTessellation/RESULT", "a+");
		fputs(line, fp);
		fclose(fp);	
    }

	delete[] voro;
}

void voronoiResult(int step)
{
	FILE * fp1;
	FILE * fp2;
	char line[512];
	int i;
	
	fp1 = fopen("voronoiTessellation/RESULT", "r");
	fp2 = fopen("voronoiTessellation/VORONOI_RESULT", "a+");	
	for(i = 0; i < NATOM; i ++){
		fgets(line, 512, fp1);
		fprintf(fp2, "Step%5d: ", step); fputs(line, fp2);
	}

	fclose(fp1);
	fclose(fp2);
}

int main(int argc, char const *argv[])
{
	FILE * fp;
	
	int i, j;
	char line[512];
	
	int start;
	start = atoi(argv[1]);

	double u[NATOM][NSTEP], v[NATOM][NSTEP], w[NATOM][NSTEP];

	mkdir("voronoiTessellation", 00777);	

	fp = fopen("XDATCAR", "r");
	for(i = 0; i < NOMIT; i ++)
		findWord("Direct configuration=", fp);
	for(i = 0; i < NSTEP; i ++){
		findWord("Direct configuration=", fp);
		fgets(line, 512, fp);
		for(j = 0; j < NATOM;j ++){
			fscanf(fp, "%lf%lf%lf", &u[j][i], &v[j][i], &w[j][i]);
		} 
	}
	fclose(fp);

	for(i = start; i < start+1000; i = i++){	

		fp = fopen("voronoiTessellation/voroimport", "w");
		for(j = 0; j < NATOM; j ++){
			fprintf(fp, "%d %lf %lf %lf\n", j+1, u[j][i], v[j][i], w[j][i]);
		}
		fclose(fp);
		
		//printf("Process step %d.\n", NOMIT+i+1);
		voronoiGenerate();
		voronoiResult(NOMIT+i+1);
	}

	if(fopen("voronoiTessellation/result1", "r") != NULL)
		remove("voronoiTessellation/result1");
	if(fopen("voronoiTessellation/result2", "r") != NULL)
		remove("voronoiTessellation/result2");
	if(fopen("voronoiTessellation/RESULT", "r") != NULL)
		remove("voronoiTessellation/RESULT");
	if(fopen("voronoiTessellation/voroimport", "r") != NULL)
		remove("voronoiTessellation/voroimport");
	
	return 0;
}
