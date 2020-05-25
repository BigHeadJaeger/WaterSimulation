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



