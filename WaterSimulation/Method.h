#pragma once
#include"GeoScene.h"
#include"mkl.h"
#include <algorithm>

//一个新的结构体用来表示映射后平面上点的信息
struct PlaneNode			//一个点映射到二维平面上后包含的信息
{
	int index;				//点的索引
	CPoint3D position;		//点在二维平面上的位置
	PlaneNode():index(0){}
};

//各种参数化方法的基类
class CMethod
{
protected:
	PlaneNode* m_Plane;
	int m_nOutPoints;					//外部点的个数
	int m_nInteralPoints;				//内部点的个数
	double** Namdar;					//系数Namdar二维数组

	vector<int> newIndex;				//前面n个是内部点的namdar，后面n+1到N是外部点的namdar,所以需要一个新的索引序列
protected:
	CMethod();
	~CMethod();

	void Initialize();					//将边界点映射到m_Plane的n+1~N位置处，并指定位置
	virtual void GetNamdar() = 0;		//根据不同的参数化方法获取Namdar
	void ConstructEquation(int flag);	//flag代表坐标点是几个数据构成的
	void Replace();						//将转换后的平面上的点替换到G_pScene
	bool judge(int v1, int v2);			//用来判断v1和v2是否邻接

	bool exist(vector<int> Array, int flag);	//判断一个数是不是在一个集合中
private:
	bool SolveEquation(vector<double> &a, vector<int> &ia, vector<int> &ja,
		vector<double> &b, vector<double> &x, int nRhs);				//解方程组（利用mkl）
};