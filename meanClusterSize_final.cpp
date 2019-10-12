#include <iostream>
#include <string.h>
#include <vector>
#include <stdio.h>

#define NSTEP 8000 
#define NATOM 200
#define Nmin 161
#define Nmax 200

using namespace std;

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

int main(int argc, char const *argv[])
{
	char line[1024];
	char neighbourstring[512];
	int id_center, temp, t, flag;
	vector<int> neighbouratom;
	vector<int> cluster_size;
	Clusters clustersPerStep;
	FILE *fp1 = fopen("voronoiTessellation/VORONOI_RESULT_total", "r");
	for(int step = 0; step < NSTEP; ++ step)
	{
		FILE *fp2 = fopen("voronoiTessellation/clusterfile", "w");
		for(int row = 0; row < NATOM; ++ row)
		{
			fgets(line, 1024, fp1);
			sscanf(line, "%*[^=]=%d%*[^=]=%*[^=]=%[^|]", &id_center, neighbourstring);
			if(id_center >= Nmin && id_center <= Nmax){
				for(int character = 0; character < strlen(neighbourstring); ++ character)
				{
					flag = sscanf(neighbourstring+character, "%d%n" ,&temp, &t); character += t;
					if(flag == 1)
						neighbouratom.push_back(temp);
				}
				fprintf(fp2, "%d ", id_center);
				for(int i = 0; i < neighbouratom.size(); ++ i)
				{
					if(neighbouratom[i] >= Nmin && neighbouratom[i] <= Nmax)
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
	cout << ratio << endl;

	FILE *fp3 = fopen("voronoiTessellation/MeanClusterSize_information", "w");
	fprintf(fp3, "%lf", ratio);
	fclose(fp3);

	return 0;
}
