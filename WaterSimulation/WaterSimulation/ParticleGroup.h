#pragma once
#include<glm.hpp>
#include<vector>
using namespace glm;
using namespace std;

#include"Object.h"
#include"PublicStruct.h"
#include"Particle.h"
#include"DistributeFun.h"
#include"MarchingCube.h"

//���ӵĽ�ģ�ӿ�(���������Ҫ���н�ģ֮������Ⱦ����metaball��)
class IModelingParticle
{
public:
	virtual void Modeling() = 0;	//�����ڵĶ��Ƿ���ֵ���ֱ�Ϊ������Ϣ���� �Ƿ��ṩ������ �Ƿ��ṩ��������
};

//����Ⱥ�Ļ���
class ParticleGroup:public Object
{
protected:
	//Render render;			//ÿһ������Ⱥ
	int particleNumber;
	vec3 centerPosition;		//����Ⱥ������λ�ã���ʱ����Ҫ
public:
};

//��ͨ����Ⱥ
class NormalParticle:public ParticleGroup
{
	
};

//MPS�㷨��ˮ����Ⱥ
class MPSWaterParticleGroup:public ParticleGroup,public IModelingParticle
{
protected:
	//�㷨������Ⱥ�������
	float l0;
	float range;
	float viscosity;
	float a;
	vector<MPSWaterParticle> particles;

	MarchingCube* marchingCube;					//һ��ָ�룬����Ҫ�õ�ʱ��ų�ʼ�������൱��һ�������
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

		marchingCube = NULL;
	}

	~MPSWaterParticleGroup()
	{
		if (marchingCube)
			delete marchingCube;
	}

	void InitParticles();
	void InitMPSTool();
	void Modeling() override;


	void Update(float dt) override;
};


