#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char const *argv[])
{
	int inputAtomID = atoi(argv[1]);
	int AtomID;
	char line[1024];
	char newline[1024];
	char string1[128], string2[128], string3[128], string4[256], string5[128];
	int temp, t, flag, length;
	vector<int> tempnumber;

	FILE *fp1; FILE *fp2;
	fp1 = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	fp2 = fopen("voronoiTessellation/Voronoi_EachAtom", "w");

	while(fgets(line, 1024, fp1))
	{
		sscanf(line, "%[^=]=%[^=]=%[^=]=%[^|]%[^\n]", string1, string2, string3, string4, string5);

		sscanf(string2, "%d%*s", &AtomID);

		if(AtomID == inputAtomID)
		{
			for(int i = 0; i < strlen(string4); i ++)
			{
				flag = sscanf(string4+i, "%d%n", &temp, &t); i+=t;
				if(flag == 1)
					tempnumber.push_back(temp);
			}
			sort(tempnumber.begin(), tempnumber.end());

			length = 0;
			for(int i = 0; i < tempnumber.size(); i ++)
			{
				sprintf(string4+length, "%d ", tempnumber[i]);
				length = strlen(string4);
			}
			sprintf(newline, "%s=%s=%s=%s%s\n", string1, string2, string3, string4 ,string5);
			fputs(newline, fp2);
			tempnumber.clear();					
		}

	}

	fclose(fp2);
	fclose(fp1);
	return 0;
}