#pragma once
#include<glm.hpp>
using namespace glm;
//����Ⱥ�Ļ���
class ParticleGroup
{
protected:
	int particleNumber;
	vec3 centerPosition;		//����Ⱥ������λ�ã���ʱ����Ҫ
public:
	virtual void Update() = 0;

};

//��ͨ����Ⱥ
class NormalParticle:public ParticleGroup
{

};

//MPS�㷨��ˮ����Ⱥ
class MPSWaterParticleGroup:public ParticleGroup
{
protected:
	//�㷨������Ⱥ�������
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