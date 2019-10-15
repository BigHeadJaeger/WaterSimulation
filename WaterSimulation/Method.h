#pragma once
#include"GeoScene.h"
#include"mkl.h"
#include <algorithm>

//һ���µĽṹ��������ʾӳ���ƽ���ϵ����Ϣ
struct PlaneNode			//һ����ӳ�䵽��άƽ���Ϻ��������Ϣ
{
	int index;				//�������
	CPoint3D position;		//���ڶ�άƽ���ϵ�λ��
	PlaneNode():index(0){}
};

//���ֲ����������Ļ���
class CMethod
{
protected:
	PlaneNode* m_Plane;
	int m_nOutPoints;					//�ⲿ��ĸ���
	int m_nInteralPoints;				//�ڲ���ĸ���
	double** Namdar;					//ϵ��Namdar��ά����

	vector<int> newIndex;				//ǰ��n�����ڲ����namdar������n+1��N���ⲿ���namdar,������Ҫһ���µ���������
protected:
	CMethod();
	~CMethod();

	void Initialize();					//���߽��ӳ�䵽m_Plane��n+1~Nλ�ô�����ָ��λ��
	virtual void GetNamdar() = 0;		//���ݲ�ͬ�Ĳ�����������ȡNamdar
	void ConstructEquation(int flag);	//flag����������Ǽ������ݹ��ɵ�
	void Replace();						//��ת�����ƽ���ϵĵ��滻��G_pScene
	bool judge(int v1, int v2);			//�����ж�v1��v2�Ƿ��ڽ�

	bool exist(vector<int> Array, int flag);	//�ж�һ�����ǲ�����һ��������
private:
	bool SolveEquation(vector<double> &a, vector<int> &ia, vector<int> &ja,
		vector<double> &b, vector<double> &x, int nRhs);				//�ⷽ���飨����mkl��
};