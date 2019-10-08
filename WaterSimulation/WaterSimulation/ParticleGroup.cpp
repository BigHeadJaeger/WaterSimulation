#include "ParticleGroup.h"

void MPSWaterParticleGroup::InitParticles()
{
	for (int i = 0; i < particleNumber; i++)
	{
		//��ʼ����ǰ���Ӳ�push��ȥ
		MPSWaterParticle R;
		//�˴�ָ���������� λ��
		R.index = i;
		particles.push_back(R);
	}
	//�����ڽӵ��ϵ
	UpdateAdjoin(range);
	//����ÿһ�����ӵĳ�ʼ�ܶ�
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
	//ÿһ֡����Ҫ�㷨���̱���ÿһ������

	//1.����ÿ�����ӵ�p
	//1.1����ÿ�����ӵ��Ҷ���
	// u*��ɢ��
	// u��������˹���
	// u*�ļ���
	// n*�ļ���
	vector<float> Right;		//�洢ÿһ�����ӵ��Ҷ���
	for (int i = 0; i < particles.size(); i++)
	{
		//������ǰ�����ӵ�λ�� �ٶȼ���
		vector<vec3> neighborPos;
		vector<vec3> neighborU;
		for (int j = 0; j < particles[i].adjoinParticleIndex.size(); j++)
		{
			neighborPos.push_back(particles[particles[i].adjoinParticleIndex[j]].position);
			neighborU.push_back(particles[particles[i].adjoinParticleIndex[j]].speed);
		}
		//������ʱ�ٶ�u* ��������ʱλ��r*
		vec3 tempU = mpsTool->TempU(dt, mpsTool->ExplicitLaplacian(neighborU, neighborPos, particles[i].speed, particles[i].position, particles[i].n0), particles[i].position);
		for (int t = 0; t < neighborPos.size(); t++)
		{
			neighborPos[t] += tempU;
		}
		vec3 tempR = particles[i].position + tempU;
		Right.push_back(mpsTool->OldImplicitLaplacianRight(particles[i].n0, dt, particles[i].n0, mpsTool->DensityN(neighborPos, tempR)));
	}
	//1.2 ������ʽ��P
	//1.2.1 ���������ӽ��б�����
	for (int i = 0; i < particles.size(); i++)
	{
		vector<vec3> neighborPos;
		neighborPos.push_back(particles[i].position);
		particles[i].SurfaceAdjudge(a, mpsTool->DensityN(neighborPos, particles[i].position), g, l0);
	}
	//1.2.2 ��һ�����ɷ���



}
