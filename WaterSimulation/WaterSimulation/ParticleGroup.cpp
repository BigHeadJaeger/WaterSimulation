#include "ParticleGroup.h"

void MPSWaterParticleGroup::InitParticles()
{
	vector<vec3> positions;
	particleNumber = CubeDistribute(positions, transformation.position, vec3(0.5), 10, 10, 10);
	for (int i = 0; i < particleNumber; i++)
	{
		//初始化当前粒子并push进去
		MPSWaterParticle R;
		//此处指定粒子属性 位置
		R.position = positions[i];
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
	vector<vec3> posArray;
	for (int i = 0; i < particles.size(); i++)
	{
		posArray.push_back(particles[i].position);
	}
	for (int i = 0; i < particles.size(); i++)
	{
		//for (int j = 0; j < particles[i].adjoinParticleIndex.size(); j++)
		//{
		//	neighborPos.push_back(particles[particles[i].adjoinParticleIndex[j]].position);
		//}
		particles[i].n0 = tool->DensityN(posArray, i);
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
	vector<double> Right;		//存储每一个粒子的右端项
	
	vector<vec3> posArray;
	vector<vec3> uArray;
	for (int i = 0; i < particles.size(); i++)
	{
		posArray.push_back(particles[i].position);
		uArray.push_back(particles[i].speed);
	}
	for (int i = 0; i < particles.size(); i++)
	{
		//声明当前点的领接点位置 速度集合

		//for (int j = 0; j < particles[i].adjoinParticleIndex.size(); j++)
		//{
		//	neighborPos.push_back(particles[particles[i].adjoinParticleIndex[j]].position);
		//	neighborU.push_back(particles[particles[i].adjoinParticleIndex[j]].speed);
		//}
		//计算临时速度u* 并更新临时位置r*
		vec3 tempU = mpsTool->TempU(dt, mpsTool->ExplicitLaplacian(uArray, posArray, i, particles[i].n0), particles[i].position);
		vector<vec3> tempPos = posArray;
		for (int t = 0; t < tempPos.size(); t++)
		{
			tempPos[t] += tempU;
		}
		vec3 tempR = particles[i].position + tempU;
		Right.push_back(mpsTool->OldImplicitLaplacianRight(particles[i].n0, dt, particles[i].n0, mpsTool->DensityN(tempPos, i)));
	}
	//1.2 计算隐式的P
	vector<bool> surfaceJudgeArray;			//记录粒子的表面判断
	vector<float> n0Array;					//记录粒子的初始密度

	//1.2.1 对所有粒子进行表面检测         疑问：此时公式中的密度怎么算，如果用之前的公式算，那是用r还是r*
	for (int i = 0; i < particles.size(); i++)
	{
		particles[i].SurfaceAdjudge(a, mpsTool->DensityN(posArray, i), g, l0);
		surfaceJudgeArray.push_back(particles[i].isSurface);
		n0Array.push_back(particles[i].n0);
	}
	//1.2.2 解一个泊松方程
	//float lambda=mpsTool->Lambda(posArray,)
	vector<double> resP = mpsTool->ImplicitCalculateP(posArray, n0Array, surfaceJudgeArray，Right);


	//计算每一个粒子新的速度U
	for (int i = 0; i < particles.size(); i++)
	{
		mat3 C = mpsTool->GetMaterixC(posArray, i, particles[i].n0);
		vec3 v1 = mpsTool->ExplicitGradient(C, resP, posArray, particles[i].n0, i);			//计算每个粒子的显式梯度
		vec3 resLU = mpsTool->ExplicitLaplacian(uArray, posArray, i, particles[i].n0);
		particles[i].speed = mpsTool->CalculateU(dt, resLU, v1, particles[i].speed, mpsTool->DensityN(posArray, i));
	}

	//更新位置
	for (int i = 0; i < particles.size(); i++)
	{
		particles[i].position += particles[i].speed;
	}

}
