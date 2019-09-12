#include "ImprovedMPS.h"

MPSToolFun* MPSToolFun::mpsTool = NULL;

template<typename T>
inline T MPSToolFun::ExplicitDivergence(vector<T> phi, vector<vec3> r, int currentIndex)
{
	return T();
}

float MPSToolFun::WeightFun(float dis)
{
	return 0.0f;
}
