#pragma once
#include<glm.hpp>
#include<vector>
using namespace glm;
using namespace std;

#include"Object.h"
#include"PublicStruct.h"
#include"Particle.h"
#include"DistributeFun.h"

//粒子群的基类
class ParticleGroup:public Object
{
protected:
	//Render render;			//每一个粒子群
	int particleNumber;
	vec3 centerPosition;		//粒子群的中心位置，暂时不需要
public:
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
	vector<MPSWaterParticle> particles;
private:
	void SetInitialN0();

	void UpdateAdjoin(float range);
public:
	MPSWaterParticleGroup()
	{
		particleNumber = 10;
		l0 = 1;
		range = 0.5;
		viscosity = 1;
		a = 0.75;
	}

	void InitParticles();
	void InitMPSTool();


	void Update(float dt) override;
};