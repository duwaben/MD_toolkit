#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

using namespace std;

#define NSTEP_voro 3000
#define NATOM 200

FILE *fp;

class Atomset
{
public:
	Atomset();
	void print();
	bool operator*(Atomset &);//判断两个集合是否有交集
	Atomset operator+(Atomset &);//求两个集合的并集
	void append(int);//添加原子
	int getsize();//得到团簇大小
	vector<int> getset();//得到团簇中的原子序号
	void reset();//团簇置零
private:
	vector<int> set;
};

//--------------------------------------------------------
Atomset::Atomset()
{
}

void Atomset::print()
{
	for(int i = 0; i < set.size(); ++ i)
	{
		cout << set[i] << " ";
	}
}

bool Atomset::operator*(Atomset &S2)//判断两个集合是否有交集
{
	for(int i = 0; i < set.size(); ++ i)
	{
		for(int j = 0; j < S2.set.size(); ++ j)
		{
			if(set[i] == S2.set[j])
				return true;
		}
	}
	return false;
}

Atomset Atomset::operator+(Atomset &S2)//求两个集合的并集
{
	Atomset temp = S2;
	int sign;
	for(int i = 0; i < set.size(); ++ i)
	{
		sign = 1;
		for(int j = 0; j < S2.set.size(); ++ j)
		{
			if(set[i] == S2.set[j])
				sign = 0;
		}
		if(sign){//如果S1中存在S2中没有的元素，则添加到temp中
			temp.append(set[i]);
		}
	}
	return temp;
}
//--------------------------------------------------------

void Atomset::append(int number)
{
	set.push_back(number);
}

int Atomset::getsize()
{
	return set.size();
}

vector<int> Atomset::getset()
{
	return set;
}

void Atomset::reset()
{
	set.clear();
}

//================================================================
//================================================================

class Clusters
{
public:
	Clusters();
	void set();
	void reset();
	void print();//输出所有团簇的大小
	vector<Atomset> getclusters();
	vector<int> getsetsize();
	void incorporationCluster(int, int);//合并团簇 
private:
	vector<Atomset> my_cluster; 
	vector<int> setsize;//纪录每个团簇的大小
};

Clusters::Clusters()
{
}

void Clusters::set()
{
	FILE* fp = fopen("voronoiTessellation/clusterfile", "r");
	char line[256];
	int temp, t, flag;
	Atomset newset;
	while(fgets(line, 256, fp)){
		for(int i = 0; i < strlen(line); i ++){
			flag = sscanf(line+i, "%d%n", &temp, &t); i += t;
			if(flag == 1)
				newset.append(temp);
		}
		my_cluster.push_back(newset);
		setsize.push_back(newset.getsize());
		newset.reset();
	}
	fclose(fp);
}

void Clusters::reset()
{
	my_cluster.clear();
	setsize.clear();
}

void Clusters::print()
{
	for(int i = 0; i < setsize.size(); ++ i){
		cout << setsize[i] << endl;
	}
}

vector<Atomset> Clusters::getclusters(){
	return my_cluster;
} 

vector<int> Clusters::getsetsize(){
	return setsize;
}

void Clusters::incorporationCluster(int num1, int num2)
{
	Atomset temp;
	temp = my_cluster[num1]+my_cluster[num2];

	my_cluster.erase(my_cluster.begin()+num1);
	setsize.erase(setsize.begin()+num1);
	my_cluster.erase(my_cluster.begin()+num2-1);
	setsize.erase(setsize.begin()+num2-1);

	my_cluster.push_back(temp);
	setsize.push_back(temp.getsize());
}

//=======================================================
//=======================================================

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

//===============================================================
//===============================================================

NumPerType readNumber(const char *name)
{
	NumPerType tempnumber;

	fp = fopen(name, "r");

	int temp, t, flag;

	const char *word1 = "Direct configuration";
	char line[512];

	//read useful informations form files
	//**********************************************************
	//**********************************************************
	
	for(int i = 0; i < 7; i ++)fgets(line, 512, fp);
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

int main(int argc, char const *argv[])
{
	//name list
	vector<string> namelist;
	readName(namelist);
	
	NumPerType thisNumPerType;
	thisNumPerType = readNumber("XDATCAR");

	//-----------------------------------------------------
	int targetAtom;
	for(int i = 0; i < namelist.size(); ++ i)
	{
		if(namelist[i] == "Re")
		{
			targetAtom = i;
		}
	}

	//-----------------------------------------------------
	char line[1024];
	char neighbourstring[512];
	int id_center, temp, t, flag;
	vector<int> neighbouratom;
	vector<int> cluster_size;
	Clusters clustersPerStep;
	FILE *fp1 = fopen("voronoiTessellation/VORONOI_RESULT", "r");
	for(int step = 0; step < NSTEP_voro; ++ step)
	{
		FILE *fp2 = fopen("voronoiTessellation/clusterfile", "w");
		for(int row = 0; row < NATOM; ++ row)
		{
			fgets(line, 1024, fp1);
			sscanf(line, "%*[^=]=%d%*[^=]=%*[^=]=%[^|]", &id_center, neighbourstring);
			if(id_center >= thisNumPerType.getmin(targetAtom)+1 && id_center < thisNumPerType.getmax(targetAtom)+1){
				for(int character = 0; character < strlen(neighbourstring); ++ character)
				{
					flag = sscanf(neighbourstring+character, "%d%n" ,&temp, &t); character += t;
					if(flag == 1)
						neighbouratom.push_back(temp);
				}
				fprintf(fp2, "%d ", id_center);
				for(int i = 0; i < neighbouratom.size(); ++ i)
				{
					if(neighbouratom[i] >= thisNumPerType.getmin(targetAtom)+1 && neighbouratom[i] < thisNumPerType.getmax(targetAtom)+1)
						fprintf(fp2, "%d ", neighbouratom[i]);
				}
				fprintf(fp2, "\n");
			}
			neighbouratom.clear();
		}
		fclose(fp2);
		//--------------------------------------------------------------
		clustersPerStep.set();
		flag = 1;
		while(flag){
			flag = 0;
			for(int i = 0; i < clustersPerStep.getclusters().size()-1 && flag == 0; ++ i)
			{
				for(int j = i+1; j < clustersPerStep.getclusters().size() && flag == 0; ++ j)
				{
					if(clustersPerStep.getclusters()[i]*clustersPerStep.getclusters()[j]){
						clustersPerStep.incorporationCluster(i, j);
						flag = 1;//找到可合并团簇
					}
				}
			}
		}
		for(int i = 0; i < clustersPerStep.getsetsize().size(); ++ i)
		{
			cluster_size.push_back(clustersPerStep.getsetsize()[i]);
		}
		clustersPerStep.reset();
	}

	fclose(fp1);

	int sumnsquare, sumn;
	sumnsquare = sumn = 0;
	for(int i = 0; i < cluster_size.size(); ++ i)
	{
		sumnsquare += cluster_size[i] * cluster_size[i];//可能出现整型溢出
		sumn += cluster_size[i];
	}
	double ratio = (double)sumnsquare /sumn;
	//cout << ratio << endl;

    //--------------------------------------------------------
	double P[thisNumPerType.getmax(targetAtom)-thisNumPerType.getmin(targetAtom)];
	for(int i = 0; i < thisNumPerType.getmax(targetAtom)-thisNumPerType.getmin(targetAtom); ++ i)
	{
		P[i] = 0;
	}

	for(int i = 0; i < thisNumPerType.getmax(targetAtom)-thisNumPerType.getmin(targetAtom); ++ i)
	{
		for(int j = 0; j < cluster_size.size(); ++ j)
		{
			if(cluster_size[j] == i+1)
			{
				P[i] = P[i] + 1.0;
			}
		}
	}

	for(int i = 0; i < thisNumPerType.getmax(targetAtom)-thisNumPerType.getmin(targetAtom); ++i)
	{
		P[i] = P[i] /cluster_size.size();
	}

	FILE *fp3 = fopen("voronoiTessellation/MeanClusterSize_information", "w");
	fprintf(fp3, "meanclustersize: %lf", ratio);
	fclose(fp3);

	fp3 = fopen("voronoiTessellation/clusterDistribution", "w");
	for(int i = 0; i < thisNumPerType.getmax(targetAtom)-thisNumPerType.getmin(targetAtom); ++ i)
	{
		fprintf(fp3, "%d\t%lf\n", i+1, P[i]);
	}
	fclose(fp3);


	//remove clusterfile
	if(fopen("voronoiTessellation/clusterfile","r") != NULL)
	{
		remove("voronoiTessellation/clusterfile");
	}

	return 0;
}