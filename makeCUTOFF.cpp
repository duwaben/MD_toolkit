#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{
	int number;
	char line[256];
	number = atoi(argv[1]);

	FILE *fp;
	fp = fopen("pairCorrelationFunctions/CUTOFFPEAK_information", "w");
	for(int i = 1; i <= number; i ++)
	{
		for(int j = i; j <= number; j ++)
		{
			sprintf(line, "CUTOFF%d%d:", i, j);
			fprintf(fp, "%s\n", line);
		}
	}
	for(int i = 1; i <= number; i ++)
	{
		for(int j = i; j <= number; j ++)
		{
			sprintf(line, "PEAK%d%d:", i, j);
			fprintf(fp, "%s\n", line);
		}
	}
	fclose(fp);
	return 0;
}
