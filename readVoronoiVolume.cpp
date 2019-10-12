#include <iostream>
#include <cstring>
#include <cstdio>
#include <vector>
#include <numeric> //accumulate

using namespace std;

#define Nmin 161
#define Nmax 200
#define NATOM 200
#define NSTEP_voro 400 

double readLattice()
{
	double lattice;
	FILE *fp = fopen("XDATCAR", "r");
	char line[512];
	fgets(line, 512, fp); fgets(line, 512, fp);
	fscanf(fp, "%lf", &lattice);
	fclose(fp);

	return lattice;
}

int main(int argc, char const *argv[])
{
	vector<double> voronoivolume;
	int id;
	double volume;
	char line[1024];
	double LATTICE = readLattice();

	FILE *fp = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	for(int i = 0; i < NATOM * NSTEP_voro; ++ i)
	{
		fgets(line, 1024, fp);
		sscanf(line, "%*[^=]=%d%*[^=]=%*[^=]=%*[^=]=%lf", &id, &volume);
		if(id >= Nmin && id <= Nmax)
		{
			voronoivolume.push_back(volume);
		}		
	}

	double sum = accumulate(voronoivolume.begin(), voronoivolume.end(), 0.0);
	double mean = sum /voronoivolume.size();
	fclose(fp);

	cout << "meanvoronoivolume =" <<mean << endl;

	fp = fopen("voronoiTessellation/MeanVoronoiVolume", "w");
	fprintf(fp, "meanvoronoivolume = %lf", mean * LATTICE * LATTICE * LATTICE);
	fclose(fp);

	return 0;
}
