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

float MPSToolFun::DensityN(vector<vec3> neighborsR, vec3 currentR)
{
	float res = 0;
	for (int i = 0; i < neighborsR.size(); i++)
	{
		res += WeightFun(distance(neighborsR[i], currentR), reForDG);
	}
	return res;
}

vec3 MPSToolFun::NewPosR(vec3 nowPos, vec3 u)
{
	return nowPos + u;
}

float MPSToolFun::Lambda(vector<vec3> neighborsR, vec3 currentR)
{
	float result = 0;
	float numerator = 0;
	float denominator = 0;
	for (int i = 0; i < neighborsR.size(); i++)
	{
		numerator += pow(distance(neighborsR[i], currentR), 2) * WeightFun(distance(neighborsR[i], currentR), reForL);
		denominator += WeightFun(distance(neighborsR[i], currentR), reForL);
	}
	result = numerator / denominator;
	return result;
}
