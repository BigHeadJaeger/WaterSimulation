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
	float Lambda(vector<vec3> r, int currentIndex);
};

template<typename T>
inline T MPSToolFun::ExplicitDivergence(vector<T> phi, vector<vec3> r, int currentIndex)
{
	T res = T();
	for (int i = 0; i < phi.size(); i++)
	{
		if (i != currentIndex)
		{
			res += ((phi[i] - phi[currentIndex]) * (r[i] - r[currentIndex]) / pow(distance(r[i], r[currentIndex]), 2))
				* WeightFun(distance(r[i], r[currentIndex]));
		}
	}
	res *= (Ds / n0);
	return res;
}

template<typename T>
inline T MPSToolFun::ExplicitLaplacian(vector<T> phi, vector<vec3> r, int currentIndex)
{
	T res = T();
	for (int i = 0; i < phi.size(); i++)
	{
		if (i != currentIndex)
		{
			res += (phi[i] - phi[currentIndex]) * WeightFun(distance(r[i], r[currentIndex]));
		}
	}

	float lambda = Lambda(r, currentIndex);

	res *= (2 * Ds / n0 * lambda);
	return res;
}

template<typename T>
inline T MPSToolFun::ExplicitGradient(vector<T> phi, vector<vec3> r, int currentIndex)
{
	return T();
}
