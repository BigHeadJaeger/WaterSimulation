#pragma once
#include<vector>
#include<glm.hpp>
using namespace glm;
using namespace std;

class MPSToolFun
{
private:
	static MPSToolFun* mpsTool;

	float n0;			//��ʼ���ܶ�ֵ
	float Ds;			//ά��
	float re;			//�ڽӵ�ķ�Χ
	MPSToolFun()
	{

	}

public:
	static MPSToolFun* GetMPSTool()
	{
		if (mpsTool == NULL)
			mpsTool = new MPSToolFun();
		return mpsTool;
	}

	//void SetN0()


	//��ʽɢ��
	template<typename T>
	T ExplicitDivergence(vector<T> phi, vector<vec3> r, int currentIndex);
	//��ʽ������˹
	template<typename T>
	T ExplicitLaplacian(vector<T> phi, vector<vec3> r, int currentIndex);
	//��ʽ�ݶ�
	template<typename T>
	T ExplicitGradient(vector<T> phi, vector<vec3> r, int currentIndex);
	//����Ȩ�صķ���
	float WeightFun(float dis);
	//MPS���ܶȵļ���
	float DensityN(vector<vec3> r, int currentIndex);
};