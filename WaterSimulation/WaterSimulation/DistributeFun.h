#pragma once
#include<vector>
#include<glm.hpp>
using namespace std;
using namespace glm;

//���� ���λ�ü��� ��ʼ��λ�� 3ά���� num�� num�� num��
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
	//���棨yֵ���ֳ�ʼ��
	for (int x = 0; x < w; x++)
		for (int z = 0; z < d; z++)
			result.push_back(initPos + vec3(offset.x * x, 0, offset.z * z));

	//�������棨һ��xֵ���ֳ�ʼ��һ��xֵ����ĩ�ˣ�
	for (int y = 0; y < h; y++)
		for (int z = 0; z < d; z++)
		{
			result.push_back(initPos + vec3(0, offset.y * y, offset.z * z));
			result.push_back(initPos + vec3(offset.x * w, offset.y * y, offset.z * z));
		}

	//����ͺ��棨һ��zֵ���ֳ�ʼ��һ��zֵ����ĩ�ˣ�
	for (int x = 0; x < h; x++)
		for (int y = 0; y < d; y++)
		{
			result.push_back(initPos + vec3(offset.x * x, offset.y * y, 0));
			result.push_back(initPos + vec3(offset.x * x, offset.y * y, offset.z * d));
		}


	return w * d + 2 * h * d + 2 * w * h - 2 * d - 2 * w - 4 * h;
}