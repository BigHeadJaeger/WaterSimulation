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



