#define _CRT_SECURE_NO_WARNINGS
#pragma once

//#include<amp.h>
#include<omp.h>
#include<math.h>
#include <stdio.h>

#include"ParticleGroup.h"
#include"MPSConst.h"

//using namespace concurrency;.

//MPS�㷨��ˮ����Ⱥ
class MPSWaterParticleGroup :public ParticleGroup, public IModelingParticle
{
protected:
	//�㷨������Ⱥ�������
	float l0;									//���ӵ�ֱ��ͬʱҲ�ǿռ��Ų�����
	float a;									//�ɳ�ϵ��
	int particleWallNum;						//ǽ���ӵĸ���
	int particleDummyNum;						//dummy�����Ӹ���
	int particleTotalNum;						//�������ӵ��ܺ�

	float n0ForNumberDensity;
	float n0ForGradient;
	float n0ForLambda;
	float lambda0;


	vector<MPSWaterParticle> particles;			//�������� wall dummy��˳�����η���

	MarchingCube* marchingCube;					//һ��ָ�룬����Ҫ�õ�ʱ��ų�ʼ�������൱��һ�������
private:
	//�õ�ǰ����Ⱥ��һЩ������ʼ��mps����
	void InitMPSTool();
	//����ÿһ�����ӵĳ�ʼ�ܶ�
	void SetInitialN0();
	////�������ӵ��ڽӹ�ϵ
	//void UpdateAdjoin(float range);
	//����Ⱥ�Ľ�ģ�ӿ�ʵ��
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

	//��ʼ������Ⱥ������
	void InitParticles();

	void ParticlesArrange();

	void Update(float dt) override;

	void Draw() override;

	void InitBufferData() override;
	void UpdateBufferData() override;
};

