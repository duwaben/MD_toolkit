#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#define NSTEP_voro 400   //8000/20
#define NATOM 200

map<string, int> dic;
map<string, int> dic1;
map<string, int> dic2;
map<string, int> dic3;
map<string, int>::iterator iter;  

void findWord(const char *word, FILE * fp)
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

int cmp(const pair<string, int>& x, const pair<string, int>& y)    
{    
    return x.second > y.second;    
}    
     
void sortMapByValue(map<string, int>& tMap,vector< pair<string, int> >& tVector)    
{    
    for (map<string, int>::iterator curr = tMap.begin(); curr != tMap.end(); curr++)     
        tVector.push_back(make_pair(curr->first, curr->second));      
     
    sort(tVector.begin(), tVector.end(), cmp);    
}    

void readNumberofAtom(const char *name, int &NATOM_1, int &NATOM_2, int &NATOM_3)
{
	FILE *fp;
	fp = fopen(name, "r");

	int i, j;

	char line[512];

	for(i = 0; i < 6; i ++)fgets(line, 512, fp);
	fscanf(fp, "%d%d%d", &NATOM_1, &NATOM_2, &NATOM_3);

	fclose(fp);
}

int main(int argc, char const *argv[])
{
	FILE *fp;
	char buffer[128];
	const char *word1 = "ID =";
	const char *word2 = "Freq_table =";
	int freq_table[7];

	//get the number of NATOM_1, NATOM_2, NATOM_3
	int NATOM_1, NATOM_2, NATOM_3;
	int atom_id, totalnumber;
	readNumberofAtom("XDATCAR", NATOM_1, NATOM_2, NATOM_3);
	//****************************************************************************************
	//****************************************************************************************
	fp = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	for(int i = 0; i < NSTEP_voro*NATOM; ++ i){
		findWord(word2, fp);
		for(int j = 0; j < 6; ++ j){
			fscanf(fp, "%d", &freq_table[j]);
		}

		if(fscanf(fp, "%d", &freq_table[6]) == 0)freq_table[6] = 0;

		sprintf(buffer, "%d %d %d %d", freq_table[3],freq_table[4],freq_table[5],freq_table[6]);
		if(dic.find(buffer) == dic.end())
			dic[buffer] = 1;
		else 
			dic[buffer] ++;
	}
	fclose(fp);

	ofstream outfile("voronoiTessellation/VoronoiDistribution_information", ios::out);
	for(iter = dic.begin(); iter != dic.end(); iter ++)  
		outfile << iter->first << "||" << iter->second << endl;
	outfile.close();

	//sort
	vector< pair<string,int> > tVector;    
    sortMapByValue(dic, tVector);
	
	fp = fopen("voronoiTessellation/VoronoiDistribution_sort_information", "w");
	for(int i = 0; i < tVector.size(); i ++)  
		fprintf(fp, "%s || %lf%%\n", tVector[i].first.c_str(), double(tVector[i].second)*100/(NSTEP_voro*NATOM));
	fclose(fp);

	//****************************************************************************************
	fp = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	totalnumber = 0;
	atom_id = 0;
	for(int i = 0; i < NSTEP_voro*NATOM; ++ i){
		
		findWord(word1, fp);
		fscanf(fp, "%d", &atom_id);
		if(atom_id > 0 && atom_id <= NATOM_1){
			totalnumber ++;
			findWord(word2, fp);
			for(int j = 0; j < 6; ++ j){
				fscanf(fp, "%d", &freq_table[j]);
			}

			if(fscanf(fp, "%d", &freq_table[6]) == 0)freq_table[6] = 0;

			sprintf(buffer, "%d %d %d %d", freq_table[3],freq_table[4],freq_table[5],freq_table[6]);
			if(dic1.find(buffer) == dic1.end())
				dic1[buffer] = 1;
			else 
				dic1[buffer] ++;
		}
	}
	fclose(fp);

	ofstream outfile1("voronoiTessellation/VoronoiDistribution_1_information", ios::out);
	for(iter = dic1.begin(); iter != dic1.end(); iter ++)  
		outfile1 << iter->first << "||" << iter->second << endl;
	outfile1.close();

	//sort
	vector< pair<string,int> > tVector1;    
    sortMapByValue(dic1, tVector1);
	
	fp = fopen("voronoiTessellation/VoronoiDistribution_1_sort_information", "w");
	for(int i = 0; i < tVector1.size(); i ++)  
		fprintf(fp, "%s || %lf%%\n", tVector1[i].first.c_str(), double(tVector1[i].second)*100/totalnumber);
	fclose(fp);

	//****************************************************************************************
	fp = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	totalnumber = 0;
	atom_id = 0;
	for(int i = 0; i < NSTEP_voro*NATOM; ++ i){
		
		findWord(word1, fp);
		fscanf(fp, "%d", &atom_id);
		if(atom_id > NATOM_1 && atom_id <= NATOM_1+NATOM_2){
			totalnumber ++;
			findWord(word2, fp);
			for(int j = 0; j < 6; ++ j){
				fscanf(fp, "%d", &freq_table[j]);
			}

			if(fscanf(fp, "%d", &freq_table[6]) == 0)freq_table[6] = 0;

			sprintf(buffer, "%d %d %d %d", freq_table[3],freq_table[4],freq_table[5],freq_table[6]);
			if(dic2.find(buffer) == dic2.end())
				dic2[buffer] = 1;
			else 
				dic2[buffer] ++;
		}
	}
	fclose(fp);

	ofstream outfile2("voronoiTessellation/VoronoiDistribution_2_information", ios::out);
	for(iter = dic2.begin(); iter != dic2.end(); iter ++)  
		outfile2 << iter->first << "||" << iter->second << endl;
	outfile2.close();

	//sort
	vector< pair<string,int> > tVector2;    
    sortMapByValue(dic2, tVector2);
	
	fp = fopen("voronoiTessellation/VoronoiDistribution_2_sort_information", "w");
	for(int i = 0; i < tVector2.size(); i ++)  
		fprintf(fp, "%s || %lf%%\n", tVector2[i].first.c_str(), double(tVector2[i].second)*100/totalnumber);
	fclose(fp);

	//****************************************************************************************
	fp = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	totalnumber = 0;
	atom_id = 0;
	for(int i = 0; i < NSTEP_voro*NATOM; ++ i){
		
		findWord(word1, fp);
		fscanf(fp, "%d", &atom_id);
		if(atom_id > NATOM_1+NATOM_2 && atom_id <= NATOM){
			totalnumber ++;
			findWord(word2, fp);
			for(int j = 0; j < 6; ++ j){
				fscanf(fp, "%d", &freq_table[j]);
			}

			if(fscanf(fp, "%d", &freq_table[6]) == 0)freq_table[6] = 0;

			sprintf(buffer, "%d %d %d %d", freq_table[3],freq_table[4],freq_table[5],freq_table[6]);
			if(dic3.find(buffer) == dic3.end())
				dic3[buffer] = 1;
			else 
				dic3[buffer] ++;
		}
	}
	fclose(fp);

	ofstream outfile3("voronoiTessellation/VoronoiDistribution_3_information", ios::out);
	for(iter = dic3.begin(); iter != dic3.end(); iter ++)  
		outfile3 << iter->first << "||" << iter->second << endl;
	outfile3.close();

	//sort
	vector< pair<string,int> > tVector3;    
    sortMapByValue(dic3, tVector3);
	
	fp = fopen("voronoiTessellation/VoronoiDistribution_3_sort_information", "w");
	for(int i = 0; i < tVector3.size(); i ++)  
		fprintf(fp, "%s || %lf%%\n", tVector3[i].first.c_str(), double(tVector3[i].second)*100/totalnumber);
	fclose(fp);
	//****************************************************************************************

	return 0;

}
