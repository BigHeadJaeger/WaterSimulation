#pragma once
#include"MPSFun.h"

//只包含最基本属性的抽象粒子
class Particle
{
public:
	vec3 position;
	bool state;
	float life;									//生命时长
	float size;									//球形时指半径
};

class MPSWaterParticle:public Particle
{
public:
	vec3 speed;									//粒子的速度向量
	float pressure;								//受到的压力
	float n0;									//初始的密度值
	bool isSurface;								//粒子是否为表面

	int index;									//记录自身的索引值
	vector<int> adjoinParticleIndex;			//记录邻接粒子的索引

public:
	MPSWaterParticle()
	{
		position = vec3(0);
		speed = vec3(0);
		pressure = 0;
		life = 0;
		size = 0;
		index = -1;
		isSurface = false;
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

	//void SetInitialN0(vector<vec3> r, int currentIndex);

	//void UpdateAdjoin(vector<Particle>& particles, float range);

	void SurfaceAdjudge(float a, float tho, float g, float l0);
};