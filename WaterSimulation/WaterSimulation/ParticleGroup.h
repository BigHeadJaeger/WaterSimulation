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
	float l0;									//���ӵ�ֱ��
	float range;								//�ڽ����ӵ��ж�����
	float viscosity;							//Һ���ճ��ϵ��
	float a;									//�ɳ�ϵ��
	vector<MPSWaterParticle> particles;

	MarchingCube* marchingCube;					//һ��ָ�룬����Ҫ�õ�ʱ��ų�ʼ�������൱��һ�������
private:
	//�õ�ǰ����Ⱥ��һЩ������ʼ��mps����
	void InitMPSTool();
	//����ÿһ�����ӵĳ�ʼ�ܶ�
	void SetInitialN0();
	//�������ӵ��ڽӹ�ϵ
	void UpdateAdjoin(float range);
	//����Ⱥ�Ľ�ģ�ӿ�ʵ��
	void Modeling() override;
public:
	MPSWaterParticleGroup()
	{
		particleNumber = 10;
		l0 = 0.01;
		range = 0.5;
		viscosity = 0.000001;
		a = 0.75;

		marchingCube = NULL;
	}

	~MPSWaterParticleGroup()
	{
		if (marchingCube)
			delete marchingCube;
	}

	//�������ӵ�ֱ��
	void SetDiameter(float d)
	{
		l0 = d;
		//ͨ��ֱ��ȷ������֮���Ӱ�췶Χ
		range = 2.1 * l0;
	}
	//����ճ��ϵ��
	void SetViscosity(float v)
	{
		viscosity = v;
	}

	//��ʼ������Ⱥ������
	void InitParticles();

	void Update(float dt) override;

	void Draw() override;

	void InitBufferData() override;
};


