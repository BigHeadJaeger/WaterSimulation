#pragma once
#include<glm.hpp>
using namespace glm;
//粒子群的基类
class ParticleGroup
{
protected:
	int particleNumber;
	vec3 centerPosition;		//粒子群的中心位置，暂时不需要
public:
	virtual void Update() = 0;

};

//普通粒子群
class NormalParticle:public ParticleGroup
{

};

//MPS算法的水粒子群
class MPSWaterParticleGroup:public ParticleGroup
{
protected:
	//算法中粒子群体的属性
	float l0;
	float range;
	float viscosity;
	float a;
public:
	MPSWaterParticleGroup()
	{
		particleNumber = 10;
		l0 = 1;
		range = 0.5;
		viscosity = 1;
		a = 0.75;
	}

	void Update() override
	{

	}
};