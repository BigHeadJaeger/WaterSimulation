//主算法进程的伪代码
#include"Particle.h"
vector<Particle> particles;

float particleNumber = 10;
float l0 = 1;
float range = 0.5;
float viscosity = 1;

int main()
{
	/*************初始化阶段****************/
	//对Mps方程工具的一些数值进行初始化
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	mpsTool->SetRe(l0);
	mpsTool->SetViscosity(viscosity);
	//粒子的初始化，各种属性的初始化，加入到粒子列表中
	for (int i = 0; i < particleNumber; i++)
	{
		//初始化当前粒子并push进去
		Particle R;
		R.index = i;
		particles.push_back(R);
	}
	//所有粒子自身的属性初始化完之后，再初始化互相之间有影响的属性
	for (int i = 0; i < particles.size(); i++)
	{
		//第一、邻接粒子的索引集合初始化
		particles[i].UpdateAdjoin(particles, range);
		//获取邻接点的位置集合
		vector<vec3> neighborPos;
		neighborPos.push_back(particles[i].position);
		for (int j = 0; j < particles[i].adjoinParticleIndex.size(); j++)
		{
			int neighborIndex = particles[i].adjoinParticleIndex[j];
			neighborPos.push_back(particles[neighborIndex].position);
		}
		//设置当前点的初始密度
		particles[i].SetInitialN0(neighborPos, particles[i].index);
	}


	/*************每一帧的流程****************/
	float deltaT = 1 / 60;
	//每一帧的主要算法流程遍历每一个粒子
	for (int i = 0; i < particles.size(); i++)
	{
		
	}



}