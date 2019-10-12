#include <iostream>
#include <cstring>
#include <cstdio>
#include <vector>
#include <numeric> //accumulate

using namespace std;

#define NATOM 200
#define NSTEP_voro 3000 

double LATTICE;

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

//-------------------------------------------------
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

void readName(vector<string> & namelist)
{
	FILE *fp;
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

NumPerType readNumber(const char *name)
{
	NumPerType tempnumber;

	FILE *fp = fopen(name, "r");

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

	fclose(fp);
	//*********************************************************
	//*********************************************************
	return tempnumber;
}

int main(int argc, char const *argv[])
{
	//namelist
	vector<string> namelist;
	readName(namelist);

	NumPerType thisNumPerType;
	thisNumPerType = readNumber("XDATCAR");
	
	vector<double> voronoivolume;
	int id;
	double volume;
	char line[1024];
	double LATTICE = readLattice();

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

	FILE *fp = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	for(int i = 0; i < NATOM * NSTEP_voro; ++ i)
	{
		fgets(line, 1024, fp);
		sscanf(line, "%*[^=]=%d%*[^=]=%*[^=]=%*[^=]=%lf", &id, &volume);
		if(id > thisNumPerType.getmin(site) && id <= thisNumPerType.getmax(site))//real serial number
		{
			voronoivolume.push_back(volume);
		}		
	}

	double sum = accumulate(voronoivolume.begin(), voronoivolume.end(), 0.0);
	double mean = sum /voronoivolume.size();
	fclose(fp);

	cout << "meanvoronoivolume =" <<mean << endl;

	fp = fopen("voronoiTessellation/MeanVoronoiVolume", "w");
	fprintf(fp, "meanvoronoivolume(relative) = %lf\n", mean);
	fprintf(fp, "meanvoronoivolume(absolute) = %lf", mean * LATTICE * LATTICE * LATTICE);
	fclose(fp);

	return 0;
}
