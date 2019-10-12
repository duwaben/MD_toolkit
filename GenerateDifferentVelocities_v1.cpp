#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <ctime>
#include <sys/stat.h>

using namespace std;

#define MASS 1.667e-27  //kg
#define kB 1.380e-23
#define Pi 3.141593
#define NATOM 200
#define NRUN 100
#define TRANS 1e5    // from m/s --> angstrom/fs

FILE *fp;
int TEBEG;

class NumPerType
{
public:
	NumPerType();
	int getmin(int);	
	int getmax(int);
	int total();
	int getsize();
	string getName(int);
private:
	vector<int> number;
	vector<string> namelist;
};

NumPerType::NumPerType()
{
	fp = fopen("POSCAR", "r");
	char line[256];
	char name[256];
	int temp, t, flag;
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

	fgets(line, 256, fp);	
	for(int i = 0; i < strlen(line); ++ i)
	{
		flag = sscanf(line+i, "%d%n", &temp, &t); i += t;
		if(flag == 1)
			number.push_back(temp);
	}

	fclose(fp);
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

int NumPerType::getsize()
{
	return number.size();
}

string NumPerType::getName(int type)
{
	return namelist[type];
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

void readINCAR(const char *name)
{
	const char *word2 = "=";
	const char *word3 = "TEBEG";

	fp = fopen(name, "r");
	findWord(word3);
	findWord(word2);
	fscanf(fp, "%d", &TEBEG);
	fclose(fp);
}

double getRMASS(const char *elementName)
{
	const char *word1 = elementName;
	const char *word2 = ":";
	fp = fopen("/public/home/users/sjtu003a/bin/RMASS_list", "r");
	findWord(word1);
	findWord(word2);
	double RMASSperType;
	fscanf(fp, "%lf", &RMASSperType);
	fclose(fp);

	return RMASSperType;
}

int main(int argc, char const *argv[])
{

	//first step: generate the intial velocities of one run.
	//second step: check if the intial velocities of two individual runs are different Yes
	//third step: large-scale manipulation
	//POSCAR INCAR RMASS_list

	srand((int)time(0));

	//generate many sets of Guassian random number
	double randn1, randn2;
	double GuassianRandomNumber[3*NATOM*NRUN];
	int numberSequence = 0;
	while(numberSequence < 3*NATOM*NRUN)
	{
		randn1 = rand();
		randn1 /= RAND_MAX;
		randn2 = rand();
		randn2 /= RAND_MAX;
		GuassianRandomNumber[numberSequence++] = sqrt(-2*log(randn2)) * cos(2*Pi*randn1); 
		GuassianRandomNumber[numberSequence++] = sqrt(-2*log(randn2)) * sin(2*Pi*randn1);
	}

	//-------------------------
	NumPerType thisNumPerType;
	readINCAR("INCAR");
	
	for(int i = 0; i < thisNumPerType.getsize(); ++i)
	{
		cout << thisNumPerType.getName(i) << "\t" << thisNumPerType.getmax(i)-thisNumPerType.getmin(i) << endl;
	}

	//arrange the intial velocities
	char filename[64];
	double RMASSperType, sigma;
	double velocity_x, velocity_y, velocity_z;
	FILE *fp1;
	numberSequence = 0;

	for(int runnumber = 1; runnumber <= 100; ++ runnumber)
	{
		sprintf(filename, "%d/InitialVelocities", runnumber);
		fp1 = fopen(filename, "w");
		for(int atomType = 0; atomType < thisNumPerType.getsize(); ++ atomType)
		{
			RMASSperType = getRMASS(thisNumPerType.getName(atomType).c_str());
			sigma = sqrt(kB*TEBEG/(RMASSperType*MASS));
			for(int j = thisNumPerType.getmin(atomType); j < thisNumPerType.getmax(atomType); ++ j)
			{
				velocity_x = sigma * GuassianRandomNumber[numberSequence++] /TRANS;
				velocity_y = sigma * GuassianRandomNumber[numberSequence++] /TRANS;
				velocity_z = sigma * GuassianRandomNumber[numberSequence++] /TRANS;
				fprintf(fp1, "%lf\t%lf\t%lf\n", velocity_x, velocity_y, velocity_z);
			}
		}

		fclose(fp1); 
	}
	
	return 0;
}