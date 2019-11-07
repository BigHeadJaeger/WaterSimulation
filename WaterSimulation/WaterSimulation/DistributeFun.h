#pragma once
#include<vector>
#include<glm.hpp>
using namespace std;
using namespace glm;

//参数 结果位置集合 初始点位置 3维步长 num宽 num高 num深
static int CubeDistribute(vector<vec3>& result, vec3 initPos, vec3 offset, int w, int h, int d)
{
	for (int x = 0; x < w; x++)
		for (int y = 0; y < h; y++)
			for (int z = 0; z < d; z++)
				result.push_back(initPos + vec3(offset.x * x, offset.y * y, offset.z * z));
	return w * h * d;
}

static int UncoverHollowCubeDistribute(vector<vec3>& result, vec3 initPos, vec3 offset, int w, int h, int d)
{
	//底面（y值保持初始）
	for (int x = 0; x < w; x++)
		for (int z = 0; z < d; z++)
			result.push_back(initPos + vec3(offset.x * x, 0, offset.z * z));

	//两个侧面（一个x值保持初始，一个x值保持末端）
	for (int y = 0; y < h; y++)
		for (int z = 0; z < d; z++)
		{
			result.push_back(initPos + vec3(0, offset.y * y, offset.z * z));
			result.push_back(initPos + vec3(offset.x * w, offset.y * y, offset.z * z));
		}

	//正面和后面（一个z值保持初始，一个z值保持末端）
	for (int x = 0; x < h; x++)
		for (int y = 0; y < d; y++)
		{
			result.push_back(initPos + vec3(offset.x * x, offset.y * y, 0));
			result.push_back(initPos + vec3(offset.x * x, offset.y * y, offset.z * d));
		}


	return w * d + 2 * h * d + 2 * w * h - 2 * d - 2 * w - 4 * h;
}