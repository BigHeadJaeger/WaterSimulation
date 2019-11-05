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

//粒子的建模接口(如果粒子需要进行建模之后再渲染，如metaball等)
class IModelingParticle
{
public:
	virtual void Modeling() = 0;	//参数内的都是返回值，分别为顶点信息数组 是否提供法向量 是否提供纹理坐标
};

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
class MPSWaterParticleGroup:public ParticleGroup,public IModelingParticle
{
protected:
	//算法中粒子群体的属性
	float l0;									//粒子的直径
	float range;								//邻接粒子的判断依据
	float viscosity;							//液体的粘度系数
	float a;									//松弛系数
	vector<MPSWaterParticle> particles;

	MarchingCube* marchingCube;					//一个指针，当需要用的时候才初始化它（相当于一个组件）
private:
	//用当前粒子群的一些参数初始化mps工具
	void InitMPSTool();
	//设置每一个粒子的初始密度
	void SetInitialN0();
	//更新粒子的邻接关系
	void UpdateAdjoin(float range);
	//粒子群的建模接口实现
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

	//设置粒子的直径
	void SetDiameter(float d)
	{
		l0 = d;
		//通过直径确定粒子之间的影响范围
		range = 2.1 * l0;
	}
	//设置粘度系数
	void SetViscosity(float v)
	{
		viscosity = v;
	}

	//初始化粒子群的属性
	void InitParticles();

	void Update(float dt) override;

	void Draw() override;

	void InitBufferData() override;
};


