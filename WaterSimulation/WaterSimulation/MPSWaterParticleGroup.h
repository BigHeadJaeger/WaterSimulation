#define _CRT_SECURE_NO_WARNINGS
#pragma once

//#include<amp.h>
#include<omp.h>
#include<math.h>
#include <stdio.h>

#include"ParticleGroup.h"
#include"MPSConst.h"

//using namespace concurrency;.

//MPS算法的水粒子群
class MPSWaterParticleGroup :public ParticleGroup, public IModelingParticle
{
protected:
	//算法中粒子群体的属性
	float l0;									//粒子的直径同时也是空间排布依据
	float a;									//松弛系数
	int particleWallNum;						//墙粒子的个数
	int particleDummyNum;						//dummy的粒子个数
	int particleTotalNum;						//所有粒子的总和

	float n0ForNumberDensity;
	float n0ForGradient;
	float n0ForLambda;
	float lambda0;


	vector<MPSWaterParticle> particles;			//按照粒子 wall dummy的顺序依次放入

	MarchingCube* marchingCube;					//一个指针，当需要用的时候才初始化它（相当于一个组件）
private:
	//用当前粒子群的一些参数初始化mps工具
	void InitMPSTool();
	//设置每一个粒子的初始密度
	void SetInitialN0();
	////更新粒子的邻接关系
	//void UpdateAdjoin(float range);
	//粒子群的建模接口实现
	void Modeling() override;





public:
	MPSWaterParticleGroup()
	{
		particleNumber = 0;
		l0 = PARTICLE_DISTANCE;
		a = 0.75;

		marchingCube = NULL;
	}

	~MPSWaterParticleGroup()
	{
		if (marchingCube)
			delete marchingCube;
	}

	//初始化粒子群的属性
	void InitParticles();

	void ParticlesArrange();

	void Update(float dt) override;

	void Draw() override;

	void InitBufferData() override;
	void UpdateBufferData() override;
};

