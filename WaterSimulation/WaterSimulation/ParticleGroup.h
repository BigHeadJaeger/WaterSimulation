#pragma once
#include<glm.hpp>
#include<vector>
using namespace glm;
using namespace std;

#include"PublicStruct.h"
#include"Particle.h"

//粒子群的基类
class ParticleGroup
{
protected:
	RenderType renderType;		//当前粒子群的渲染方式
	//Render render;			//每一个粒子群
	int particleNumber;
	vec3 centerPosition;		//粒子群的中心位置，暂时不需要
public:
	virtual void Update() = 0;	//抽象更新粒子状态的函数，不同的粒子群有不同的更新方式

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
		renderType = PARTICLE_BILL_BOARD;
		particleNumber = 10;
		l0 = 1;
		range = 0.5;
		viscosity = 1;
		a = 0.75;
	}

	void InitParticles();
	void InitMPSTool();


	void Update() override
	{

	}
};