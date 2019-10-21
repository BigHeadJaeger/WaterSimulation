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
	float Lambda(vector<vec3> r, int currentIndex);
	//����Ȩ�صķ���
	float WeightFun(float dis, float re);
	//������ʱֵu* (ʱ��u��������˹�������ǰ���ӵ�u)
	vec3 TempU(float deltaT, vec3 resLU, vec3 uNow);
	//������ʽ������˹���Ҷ���(��ʼ�ܶȣ�u*��ɢ�ȣ�ʱ��n0��n*)
	vec3 ImplicitLaplacianRight(float rho0, vec3 resDu, float deltaT, float n0, float tempN);
	float OldImplicitLaplacianRight(float rho0, float deltaT, float n0, float tempN);
	//��ȡ����C
	mat3 GetMaterixC(vector<vec3>& R, int currentIndex, float n0);

	//�ʹ��ͷ�������صĺ���
	//1.���̵Ĺ���
	//void ConstructEquation(int dimension);

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
	T ExplicitDivergence(vector<T>& phi, vector<vec3>& r, int currentIndex, float n0);
	//��ʽ������˹
	template<typename T>
	T ExplicitLaplacian(vector<T>& phi, vector<vec3>& r, int currentIndex, float n0);
	//��ʽ�ݶȣ��·��������ݶȣ�
	vec3 ExplicitGradient(mat3 C, vector<double>& p, vector<vec3>& r, float n0, int currentIndex);

	//MPS���ܶȵļ���
	float DensityN(vector<vec3> r, int currentIndex);
	//�����ٶ���u�����µ�λ��
	vec3 NewPosR(vec3 nowPos, vec3 u);
	//������ʵ��uֵ
	vec3 CalculateU(float deltaT, vec3 resLU, vec3 resGP, vec3 uNow, float tho);
	//��ʽ����P����ϡ�跽���飩
	vector<double> ImplicitCalculateP(vector<vec3>& r, vector<float>& n0, vector<bool>& isSurface);
};

template<typename T>
inline T MPSToolFun::ExplicitDivergence(vector<T>& phi, vector<vec3>& r, int currentIndex, float n0)
{
	T res = T();
	for (int i = 0; i < phi.size(); i++)
	{
		if (i != currentIndex)
		{
			res += ((phi[i] - phi[currentIndex]) * (r[i] - r[currentIndex]) / pow(distance(r[i], r[currentIndex]), 2))
				* WeightFun(distance(r[i], r[currentIndex]), reForDG);
		}
	}
	res *= (Ds / n0);
	return res;
}

template<typename T>
inline T MPSToolFun::ExplicitLaplacian(vector<T>& phi, vector<vec3>& r, int currentIndex, float n0)
{
	T res = T();
	for (int i = 0; i < phi.size(); i++)
	{
		if (i != currentIndex)
		{
			res += (phi[i] - phi[currentIndex]) * WeightFun(distance(r[i], r[currentIndex]), reForL);
		}
	}

	float lambda = Lambda(r, currentIndex);

	res *= (2 * Ds / n0 * lambda);
	return res;
}
