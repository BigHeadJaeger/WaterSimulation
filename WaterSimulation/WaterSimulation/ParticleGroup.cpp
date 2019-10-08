#include "ParticleGroup.h"

void MPSWaterParticleGroup::InitParticles()
{
	for (int i = 0; i < particleNumber; i++)
	{
		//初始化当前粒子并push进去
		MPSWaterParticle R;
		//此处指定粒子属性 位置
		R.index = i;
		particles.push_back(R);
	}
	//更新邻接点关系
	UpdateAdjoin(range);
	//设置每一个粒子的初始密度
	SetInitialN0();


}

void MPSWaterParticleGroup::InitMPSTool()
{
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	mpsTool->SetRe(l0);
	mpsTool->SetViscosity(viscosity);
}

void MPSWaterParticleGroup::SetInitialN0()
{
	MPSToolFun* tool = MPSToolFun::GetMPSTool();
	for (int i = 0; i < particles.size(); i++)
	{
		vector<vec3> neighborPos;
		for (int j = 0; j < particles[i].adjoinParticleIndex.size(); j++)
		{
			neighborPos.push_back(particles[particles[i].adjoinParticleIndex[j]].position);
		}
		particles[i].n0 = tool->DensityN(neighborPos, particles[i].position);
	}

}

void MPSWaterParticleGroup::UpdateAdjoin(float range)
{
	for (int i = 0; i < particles.size(); i++)
	{
		particles[i].adjoinParticleIndex.clear();
		for (int j = 0; j < particles.size(); j++)
		{
			if (j != i)
			{
				vec3 pos1 = particles[j].position;
				if (distance(pos1, particles[i].position) <= range)
					particles[i].adjoinParticleIndex.push_back(j);
			}
		}
	}

}

void MPSWaterParticleGroup::Update(float dt)
{
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	//每一帧的主要算法流程遍历每一个粒子

	//1.计算每个粒子的p
	//1.1计算每个粒子的右端项
	// u*的散度
	// u的拉普拉斯结果
	// u*的计算
	// n*的计算
	vector<float> Right;		//存储每一个粒子的右端项
	for (int i = 0; i < particles.size(); i++)
	{
		//声明当前点的领接点位置 速度集合
		vector<vec3> neighborPos;
		vector<vec3> neighborU;
		for (int j = 0; j < particles[i].adjoinParticleIndex.size(); j++)
		{
			neighborPos.push_back(particles[particles[i].adjoinParticleIndex[j]].position);
			neighborU.push_back(particles[particles[i].adjoinParticleIndex[j]].speed);
		}
		//计算临时速度u* 并更新临时位置r*
		vec3 tempU = mpsTool->TempU(dt, mpsTool->ExplicitLaplacian(neighborU, neighborPos, particles[i].speed, particles[i].position, particles[i].n0), particles[i].position);
		for (int t = 0; t < neighborPos.size(); t++)
		{
			neighborPos[t] += tempU;
		}
		vec3 tempR = particles[i].position + tempU;
		Right.push_back(mpsTool->OldImplicitLaplacianRight(particles[i].n0, dt, particles[i].n0, mpsTool->DensityN(neighborPos, tempR)));
	}
	//1.2 计算隐式的P
	//1.2.1 对所有粒子进行表面检测
	for (int i = 0; i < particles.size(); i++)
	{
		vector<vec3> neighborPos;
		neighborPos.push_back(particles[i].position);
		particles[i].SurfaceAdjudge(a, mpsTool->DensityN(neighborPos, particles[i].position), g, l0);
	}
	//1.2.2 解一个泊松方程



}
