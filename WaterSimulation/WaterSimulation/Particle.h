#pragma once
//#include"MPSFun.h"
#include<vector>
#include<glm.hpp>
#include<iostream>
using namespace glm;
using namespace std;

//只包含最基本属性的抽象粒子
class Particle
{
public:
	vec3 position;
	bool state;
	float life;									//生命时长
	float size;									//球形时指半径
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
	vec3 speed;									//粒子的速度向量
	float pressure;								//受到的压力

	float n0;									//初始的密度值
	float middleN;
	// 三种初始密度n0
	float n0ForNumberDensity;
	float n0ForGradient;
	float n0ForLambda;


	float tho;									//当前的密度


	//ParticleType particleType;
	//bool isFluid;								//是流体
	//bool isSurface;								//粒子是否为表面
	//bool isWall;								//粒子是否为墙
	//bool isDummy;								//粒子是否为虚拟墙面（作用是区分墙粒子和表面粒子）

	int index;									//记录自身的索引值
	vector<int> adjoinParticleIndex;			//记录邻接粒子的索引

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