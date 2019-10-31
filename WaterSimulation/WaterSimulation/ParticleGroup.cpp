#include "ParticleGroup.h"

void MPSWaterParticleGroup::InitParticles()
{
	vector<vec3> positions;
	particleNumber = CubeDistribute(positions, transformation.position, vec3(0.5), 10, 10, 10);
	for (int i = 0; i < particleNumber; i++)
	{
		//��ʼ����ǰ���Ӳ�push��ȥ
		MPSWaterParticle R;
		//�˴�ָ���������� λ��
		R.position = positions[i];
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

void MPSWaterParticleGroup::Modeling()
{
	vector<vec3> posArray;
	for (int i = 0; i < particles.size(); i++)
	{
		posArray.push_back(particles[i].position);
	}

	vector<float> verticesInfo;
	bool provideNormal;
	bool provideTex;

	MarchingCube::GetInstance()->GetMeshData(posArray, verticesInfo, provideNormal, provideTex);

	shaderData->InitVertexBuffer(verticesInfo, provideNormal, provideTex);
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
	//ÿһ֡����Ҫ�㷨���̱���ÿһ������

	//1.����ÿ�����ӵ�p
	//1.1����ÿ�����ӵ��Ҷ���
	// u*��ɢ��
	// u��������˹���
	// u*�ļ���
	// n*�ļ���
	vector<double> Right;					//�洢ÿһ�����ӵ��Ҷ���
	
	vector<vec3> posArray;					//��¼kʱ���������ӵ�λ��
	vector<vec3> uArray;					//��¼kʱ���������ӵ��ٶ�
	vector<vec3> tempUArray;				//��¼k+1ʱ���������ӵ�u*
	vector<vec3> tempPosArray;				//��¼k+1ʱ���������ӵ�r*

	//1.2 ������ʽ��P
	vector<bool> surfaceJudgeArray;			//��¼���ӵı����ж�
	vector<float> n0Array;					//��¼���ӵĳ�ʼ�ܶ�

	for (int i = 0; i < particles.size(); i++)
	{
		posArray.push_back(particles[i].position);
		uArray.push_back(particles[i].speed);
	}
	//�������� ��ȡ��ʱ���ٶȺ�λ��
	for (int i = 0; i < particles.size(); i++)
	{
		//������ʱ�ٶ�u* ��������ʱλ��r*
		vec3 tempU = mpsTool->TempU(dt, mpsTool->ExplicitLaplacian(uArray, posArray, i, particles[i].n0), particles[i].position);
		tempUArray.push_back(tempU);

		vec3 tempPos = particles[i].position + tempU;
		tempPosArray.push_back(tempPos);
	}
	//����ʱ�ٶȺ���ʱλ�ü���ÿ�����Ӷ�Ӧ���Ҷ���
	for (int i = 0; i < particles.size(); i++)
	{
		float resDivergence = mpsTool->ExplicitDivergence(tempUArray, tempPosArray, i, particles[i].n0);
		float resRight = mpsTool->ImplicitLaplacianRight(particles[i].n0, resDivergence, dt, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
		Right.push_back(resRight);

		//��Ϊ���������Ҷ���ļ��㲻Ӱ�죬����ͬһ��ѭ�������Ч��
		particles[i].SurfaceAdjudge(a, mpsTool->DensityN(posArray, i), g, l0);
		surfaceJudgeArray.push_back(particles[i].isSurface);
		n0Array.push_back(particles[i].n0);
	}

	//1.2.2 ��һ�����ɷ���
	vector<double> resP = mpsTool->ImplicitCalculateP(posArray, n0Array, surfaceJudgeArray, Right);


	//����ÿһ�������µ��ٶ�U
	for (int i = 0; i < particles.size(); i++)
	{
		mat3 C = mpsTool->GetMaterixC(posArray, i, particles[i].n0);
		vec3 v1 = mpsTool->ExplicitGradient(C, resP, posArray, particles[i].n0, i);			//����ÿ�����ӵ���ʽ�ݶ�
		vec3 resLU = mpsTool->ExplicitLaplacian(uArray, posArray, i, particles[i].n0);
		particles[i].speed = mpsTool->CalculateU(dt, resLU, v1, particles[i].speed, mpsTool->DensityN(posArray, i));
	}

	//����λ��
	for (int i = 0; i < particles.size(); i++)
	{
		particles[i].position += particles[i].speed;
	}

	//������õ�λ�õ���н�ģ
	//Modeling();
}
