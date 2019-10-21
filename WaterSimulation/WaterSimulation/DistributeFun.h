#pragma once
#include<vector>
#include<glm.hpp>
using namespace std;
using namespace glm;

//参数 结果位置集合 初始点位置 3维步长 num宽 num高 num深
int CubeDistribute(vector<vec3>& result, vec3 initPos, vec3 offset, int w, int h, int d)
{
	for (int x = 0; x < w; x++)
		for (int y = 0; y < h; y++)
			for (int z = 0; z < d; z++)
				result.push_back(initPos + vec3(offset.x * x, offset.y * y, offset.z * z));
	return w * h * d;
}