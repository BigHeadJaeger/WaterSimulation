#pragma once
#include<vector>
#include<glm.hpp>
using namespace glm;
using namespace std;

#include"Const.h"

class MPSToolFun
{
private:
	static MPSToolFun* mpsTool;


	float Ds;				//ά��
	float reForDG;			//�����ܶȺ��ݶȵ�re
	float reForL;			//����������˹��re
	float viscosity;		//ճ��ϵ��
	float gama;				//
	MPSToolFun()
	{
		Ds = 3;
		reForDG = 0;
		reForL = 0;
		viscosity = 0;
		gama = 0.01;
	}

public:
	//�ڲ�˽�м��㺯��
	//Lambda����
	float Lambda(vector<vec3> neighborsR, vec3 currentR);
	//����Ȩ�صķ���
	float WeightFun(float dis, float re);
	//������ʱֵu* (ʱ��u��������˹�������ǰ���ӵ�u)
	vec3 TempU(float deltaT, vec3 resLU, vec3 uNow);
	//������ʽ������˹���Ҷ���(��ʼ�ܶȣ�u*��ɢ�ȣ�ʱ��n0��n*)
	vec3 ImplicitLaplacianRight(float rho0, vec3 resDu, float deltaT, float n0, float tempN);
	float OldImplicitLaplacianRight(float rho0, float deltaT, float n0, float tempN);
public:
	//�ⲿ�ӿ�
	static MPSToolFun* GetMPSTool()
	{
		if (mpsTool == NULL)
			mpsTool = new MPSToolFun();
		return mpsTool;
	}

	//void SetN0()
	void SetRe(float l0)
	{
		reForDG = 2.1 * l0;
		reForL = 3.1 * l0;
	}

	void SetViscosity(float v)
	{
		viscosity = v;
	}

	//��ʽɢ��
	template<typename T>
	T ExplicitDivergence(vector<T>& neighborsPhi, vector<vec3>& neighborsR, T nowPhi, vec3 nowR, float n0);
	//��ʽ������˹
	template<typename T>
	T ExplicitLaplacian(vector<T>& neighborsPhi, vector<vec3>& neighborsR, T nowPhi, vec3 nowR, float n0);
	//��ʽ�ݶ�
	template<typename T>
	T ExplicitGradient(vector<T> phi, vector<vec3> r, int currentIndex);

	//MPS���ܶȵļ���
	float DensityN(vector<vec3> neighborsR, vec3 currentR);
	//�����ٶ���u�����µ�λ��
	vec3 NewPosR(vec3 nowPos, vec3 u);
};

template<typename T>
inline T MPSToolFun::ExplicitDivergence(vector<T>& neighborsPhi, vector<vec3>& neighborsR, T nowPhi, vec3 nowR, float n0)
{
	T res = T();
	for (int i = 0; i < neighborsPhi.size(); i++)
	{
		res += ((neighborsPhi[i] - nowPhi) * (neighborsR[i] - nowR) / pow(distance(neighborsR[i], nowR), 2))
			* WeightFun(distance(neighborsR[i], nowR), reForDG);
	}
	res *= (Ds / n0);
	return res;
}

template<typename T>
inline T MPSToolFun::ExplicitLaplacian(vector<T>& neighborsPhi, vector<vec3>& neighborsR, T nowPhi, vec3 nowR, float n0)
{
	T res = T();
	for (int i = 0; i < neighborsPhi.size(); i++)
	{
		res += (neighborsPhi[i] - nowPhi) * WeightFun(distance(neighborsR[i], nowR), reForL);
	}

	float lambda = Lambda(neighborsR, nowR);

	res *= (2 * Ds / n0 * lambda);
	return res;
}

template<typename T>
inline T MPSToolFun::ExplicitGradient(vector<T> phi, vector<vec3> r, int currentIndex)
{
	return T();
}
