#include <iostream>
#include <fstream>
#include <string.h>
using namespace std;

int main(int argc, char const *argv[])
{
	fstream fp, outfile;
	char line[256];
	int v1, v2, v3, v4, number, totalnumber;
	int hist[21];
	//
	
	for(int i = 0; i < 21; ++ i)
		hist[i] = 0;
	totalnumber = 0;
	fp.open("voronoiTessellation/VoronoiDistribution_information", ios::in);
	while(fp.getline(line, 256)){
		sscanf(line, "%d%d%d%d||%d", &v1, &v2, &v3, &v4, &number);		
		hist[v1+v2+v3+v4] += number;
		totalnumber += number;
	}
	fp.close(); 	
	
	outfile.open("voronoiTessellation/CNdistribution_information", ios::out);
	for(int i = 0; i < 21; ++ i)
		outfile << i << "\t" << double(hist[i])/totalnumber*100 << "%" <<endl;
	outfile.close();

	for(int i = 0; i < 21; ++ i)
		hist[i] = 0;
	totalnumber = 0;
	fp.open("voronoiTessellation/VoronoiDistribution_1_information", ios::in);
	while(fp.getline(line, 256)){
		sscanf(line, "%d%d%d%d||%d", &v1, &v2, &v3, &v4, &number);		
		hist[v1+v2+v3+v4] += number;
		totalnumber += number;
	}
	fp.close(); 	
	
	outfile.open("voronoiTessellation/CNdistribution_1_information", ios::out);
	for(int i = 0; i < 21; ++ i)
		outfile << i << "\t" << double(hist[i])/totalnumber*100 << "%" <<endl;
	outfile.close();

	for(int i = 0; i < 21; ++ i)
		hist[i] = 0;
	totalnumber = 0;
	fp.open("voronoiTessellation/VoronoiDistribution_2_information", ios::in);
	while(fp.getline(line, 256)){
		sscanf(line, "%d%d%d%d||%d", &v1, &v2, &v3, &v4, &number);		
		hist[v1+v2+v3+v4] += number;
		totalnumber += number;
	}
	fp.close(); 	
	
	outfile.open("voronoiTessellation/CNdistribution_2_information", ios::out);
	for(int i = 0; i < 21; ++ i)
		outfile << i << "\t" << double(hist[i])/totalnumber*100 << "%" <<endl;
	outfile.close();


	for(int i = 0; i < 21; ++ i)
		hist[i] = 0;
	totalnumber = 0;
	fp.open("voronoiTessellation/VoronoiDistribution_3_information", ios::in);
	while(fp.getline(line, 256)){
		sscanf(line, "%d%d%d%d||%d", &v1, &v2, &v3, &v4, &number);		
		hist[v1+v2+v3+v4] += number;
		totalnumber += number;
	}
	fp.close(); 	
	
	outfile.open("voronoiTessellation/CNdistribution_3_information", ios::out);
	for(int i = 0; i < 21; ++ i)
		outfile << i << "\t" << double(hist[i])/totalnumber*100 << "%" <<endl;
	outfile.close();


	return 0;
}
