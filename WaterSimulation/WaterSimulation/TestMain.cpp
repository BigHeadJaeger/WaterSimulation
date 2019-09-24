//主算法进程的伪代码
#include"Particle.h"
vector<Particle> particles;

float particleNumber = 10;
float l0 = 1;
float range = 0.5;
float viscosity = 1;
float a = 0.75;


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

	//1.计算每个粒子的p
	//1.1计算每个粒子的右端项
	// u*的散度
	// u的拉普拉斯结果
	// u*的计算
	// n*的计算
	vector<float> Right;		//存储每一个粒子的右端项
	for (int i = 0; i < particles.size(); i++)
	{
		vector<vec3> neighborPos;
		vector<vec3> neighborU;
		neighborPos.push_back(particles[i].position);
		neighborU.push_back(particles[i].speed);
		vec3 tempU = mpsTool->TempU(deltaT, mpsTool->ExplicitLaplacian(neighborU, neighborPos, i, particles[i].n0), particles[i].position);
		Right.push_back(mpsTool->OldImplicitLaplacianRight(particles[i].n0, deltaT, particles[i].n0, mpsTool->DensityN(neighborPos, i)));
	}
	//1.2 计算隐式的P
	//1.2.1 对所有粒子进行表面检测
	for (int i = 0; i < particles.size(); i++)
	{
		vector<vec3> neighborPos;
		neighborPos.push_back(particles[i].position);
		particles[i].SurfaceAdjudge(a, mpsTool->DensityN(neighborPos, i), g, l0);
	}
	//1.2.2 解一个泊松方程



}