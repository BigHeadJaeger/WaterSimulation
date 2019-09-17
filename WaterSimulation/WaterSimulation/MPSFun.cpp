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

float MPSToolFun::Lambda(vector<vec3> r, int currentIndex)
{
	float result = 0;
	float numerator = 0;
	float denominator = 0;
	for (int i = 0; i < r.size(); i++)
	{
		if (i != currentIndex)
		{
			numerator += pow(distance(r[i], r[currentIndex]), 2) * WeightFun(distance(r[i], r[currentIndex]));
			denominator += WeightFun(distance(r[i], r[currentIndex]));
		}
	}
	result = numerator / denominator;
	return result;
}
