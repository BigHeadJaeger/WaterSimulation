#pragma once
#include<glm.hpp>
#include<vector>
using namespace glm;
using namespace std;

#include"PublicStruct.h"
#include"Particle.h"

//����Ⱥ�Ļ���
class ParticleGroup
{
protected:
	RenderType renderType;		//��ǰ����Ⱥ����Ⱦ��ʽ
	//Render render;			//ÿһ������Ⱥ
	int particleNumber;
	vec3 centerPosition;		//����Ⱥ������λ�ã���ʱ����Ҫ
public:
	virtual void Update() = 0;	//�����������״̬�ĺ�������ͬ������Ⱥ�в�ͬ�ĸ��·�ʽ

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