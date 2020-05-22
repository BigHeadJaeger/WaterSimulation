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
#include"Const.h"

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


