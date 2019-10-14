#include "MPSFun.h"

MPSToolFun* MPSToolFun::mpsTool = NULL;


float MPSToolFun::WeightFun(float dis,float re)
{
	float result = 0;
	if (dis >= 0 && dis < re)
		result = re / dis - 1;
	else if (dis >= re)
		result = 0;
	return result;
}

vec3 MPSToolFun::TempU(float deltaT, vec3 resLU, vec3 uNow)
{
	return deltaT * (viscosity * resLU + g) + uNow;
}

vec3 MPSToolFun::ImplicitLaplacianRight(float rho0, vec3 resDu, float deltaT, float n0, float tempN)
{
	return (1 - gama) * rho0 * (resDu / deltaT) - gama * (rho0 / pow(deltaT, 2)) * ((tempN - n0) / n0);
}

float MPSToolFun::OldImplicitLaplacianRight(float rho0, float deltaT, float n0, float tempN)
{
	return (rho0 / pow(deltaT, 2)) * ((n0 - tempN) / n0);
}

mat3 MPSToolFun::GetMaterixC(vector<vec3> R, int currentIndex, float n0)
{
	mat3 res(0);
	for (int i = 0; i < R.size(); i++)
	{
		if (i != currentIndex)
		{
			vec3 v1 = (1 / length(R[i] - R[currentIndex])) * (R[i] - R[currentIndex]);
			res += (WeightFun(length(R[i] - R[currentIndex]), reForDG) * outerProduct(v1, v1));
		}
	}
	res *= (1 / n0);
	return res;
}

vec3 MPSToolFun::ExplicitGradient(mat3 C, vector<float>& p, vector<vec3>& r, float n0, float Ds,int currentIndex)
{
	if (determinant(C) >= 0.05)
	{
		vec3 res(0);
		for (int i = 0; i < p.size(); i++)
		{
			if (i != currentIndex)
			{
				float l = length(r[i] - r[currentIndex]);
				res += (WeightFun(length(l), reForDG) * ((p[i] - p[currentIndex]) / l) * ((r[i] - r[currentIndex]) / l));
			}
		}
		res /= n0;
		return inverse(C)* res;
	}
	else
	{
		vec3 res(0);
		for (int i = 0; i < p.size(); i++)
		{
			if (i != currentIndex)
			{
				float l = length(r[i] - r[currentIndex]);
				res += (WeightFun(length(l), reForDG) * ((p[i] - p[currentIndex]) / l) * ((r[i] - r[currentIndex]) / l));
			}
		}
		res *= (Ds / n0);
		return res;
	}
}

float MPSToolFun::DensityN(vector<vec3> r, int currentIndex)
{
	float res = 0;
	for (int i = 0; i < r.size(); i++)
	{
		if (i != currentIndex)
		{
			res += WeightFun(distance(r[i], r[currentIndex]), reForDG);
		}
	}
	return res;
}

vec3 MPSToolFun::NewPosR(vec3 nowPos, vec3 u)
{
	return nowPos + u;
}

float MPSToolFun::Lambda(vector<vec3> r, int currentIndex)
{
	float result = 0;
	float numerator = 0;
	float denominator = 0;
	for (int i = 0; i < r.size(); i++)
	{
		if (i != currentIndex)
		{
			numerator += pow(distance(r[i], r[currentIndex]), 2) * WeightFun(distance(r[i], r[currentIndex]), reForL);
			denominator += WeightFun(distance(r[i], r[currentIndex]), reForL);
		}
	}
	result = numerator / denominator;
	return result;
}
