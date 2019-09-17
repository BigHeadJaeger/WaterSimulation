#pragma once
#include<vector>
#include<glm.hpp>
using namespace glm;
using namespace std;

class MPSToolFun
{
private:
	static MPSToolFun* mpsTool;


	float Ds;			//维数
	float reForDG;			//计算密度和梯度的re
	float reForL;			//计算拉普拉斯的re
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
	void SetRe(float l0)
	{
		reForDG = 2.1 * l0;
		reForL = 3.1 * l0;
	}

	//显式散度
	template<typename T>
	T ExplicitDivergence(vector<T> phi, vector<vec3> r, int currentIndex);
	//显式拉普拉斯
	template<typename T>
	T ExplicitLaplacian(vector<T> phi, vector<vec3> r, int currentIndex);
	//显式梯度
	template<typename T>
	T ExplicitGradient(vector<T> phi, vector<vec3> r, int currentIndex);
	//计算权重的方程
	float WeightFun(float dis, float re);
	//MPS中密度的计算
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
				* WeightFun(distance(r[i], r[currentIndex]), reForDG);
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
			res += (phi[i] - phi[currentIndex]) * WeightFun(distance(r[i], r[currentIndex]), reForL);
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
