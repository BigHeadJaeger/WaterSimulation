#pragma once
//#include"MPSFun.h"
#include<vector>
#include<glm.hpp>
#include<iostream>
using namespace glm;
using namespace std;

//ֻ������������Եĳ�������
class Particle
{
public:
	vec3 position;
	bool state;
	float life;									//����ʱ��
	float size;									//����ʱָ�뾶
};

enum ParticleType
{
	ASDASD,
	//FLUID,
	//SURFACE,
	//WALL,
	//DUMMY
};

class MPSWaterParticle:public Particle
{
public:
	vec3 speed;									//���ӵ��ٶ�����
	float pressure;								//�ܵ���ѹ��

	float n0;									//��ʼ���ܶ�ֵ
	float middleN;
	// ���ֳ�ʼ�ܶ�n0
	float n0ForNumberDensity;
	float n0ForGradient;
	float n0ForLambda;


	float tho;									//��ǰ���ܶ�


	//ParticleType particleType;
	//bool isFluid;								//������
	//bool isSurface;								//�����Ƿ�Ϊ����
	//bool isWall;								//�����Ƿ�Ϊǽ
	//bool isDummy;								//�����Ƿ�Ϊ����ǽ�棨����������ǽ���Ӻͱ������ӣ�

	int index;									//��¼���������ֵ
	vector<int> adjoinParticleIndex;			//��¼�ڽ����ӵ�����

public:
	MPSWaterParticle()
	{
		position = vec3(0);
		speed = vec3(0);
		pressure = 0;
		n0 = 0;
		tho = 0;
		life = 0;
		size = 0;
		index = -1;
		//particleType = FLUID;
		middleN = 0;
		//isFluid = true;
		//isSurface = false;
		//isWall = false;
		//isDummy = false;
	}
	MPSWaterParticle(vec3 _pos, vec3 _speed, float _pressure, float _size = 1, float _life = 0)
	{
		position = _pos;
		speed = _speed;
		pressure = _pressure;
		life = _life;
		size = _size;
		index = -1;
	}

	//void SurfaceAdjudge(float a, float g, float l0, float thoNow);
	//void OldSurfaceAdjudge(float b, float nNow);
};