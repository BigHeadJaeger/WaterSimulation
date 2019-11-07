#include "ParticleGroup.h"

void MPSWaterParticleGroup::InitParticles()
{
	//��ʼ��MPS����
	InitMPSTool();
	vec3 initPos = transformation.position;
	vec3 offset = vec3(l0);									//���ӿռ�ֲ�����
	int width, height, depth;								//���ӵĳ�����Ų�
	width = height = depth = 10;

	float highest = 0;										//�ҵ���ߵ�λ��
	if (offset.y > 0)
		highest = initPos.y + (height - 1) * offset.y;
	else
		highest = initPos.y;
	vector<vec3> positions;
	int n = 0;
	n = CubeDistribute(positions, initPos, offset, width, height, depth);
	for (int i = 0; i < n; i++)
	{
		//��ʼ����ǰ���Ӳ�push��ȥ
		MPSWaterParticle R;
		//�˴�ָ���������� λ��
		R.position = positions[i];
		//����������������������Ϊ��������
		if (positions[i].y == highest)
			R.isSurface = true;
		R.index = i;
		particles.push_back(R);
	}
	particleNumber = n;

	//����Ҫ��ʼ������ǽ�Լ�dummy wall	(���Կ����ٴ���1+2������޸ǵ�cube)   ps:Ҫȷ��������Χס���е�����
	//wall
	int containerWidth = width * 2;
	int containerHeight = height + 10;
	int containerDepth = depth + 2;
	initPos -= offset;
	positions.clear();
	n = UncoverHollowCubeDistribute(positions, initPos, offset, containerWidth, containerHeight, containerDepth);
	for (int i = 0; i < n; i++)
	{
		MPSWaterParticle R;
		R.position = positions[i];
		R.isWall = true;
		R.index = i + particleNumber;				//����֮ǰ�Ѿ��е�����ƫ����
	}
	particleWallNum = n;
	//dummy
	int dummyWidth = containerWidth + 2;
	int dummyHeight = containerHeight;
	int dummyDepth = containerDepth + 2;
	positions.clear();

	initPos -= offset;
	n = UncoverHollowCubeDistribute(positions, initPos, offset, containerWidth, containerHeight, containerDepth);
	//�ڶ���
	dummyWidth += 2;
	dummyDepth += 2;
	initPos -= offset;
	n += UncoverHollowCubeDistribute(positions, initPos, offset, containerWidth, containerHeight, containerDepth);
	for (int i = 0; i < n; i++)
	{
		MPSWaterParticle R;
		R.position = positions[i];
		R.isDummy = true;
		R.index = i + particleNumber + particleWallNum;				//����֮ǰ�Ѿ��е�����ƫ����
	}
	particleDummyNum = n;

	particleTotalNum = particleNumber + particleDummyNum + particleWallNum;

	//�����ڽӵ��ϵ
	UpdateAdjoin(range);
	//����ÿһ�����ӵĳ�ʼ�ܶ�
	SetInitialN0();


	//�˴������������ӵĳ�ʼѹ��
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	vector<double> Right;					//�洢ÿһ�����ӵ��Ҷ���

	vector<vec3> posArray;					//��¼kʱ���������ӵ�λ��
	vector<vec3> posArray1;					//��¼kʱ�̳�ȥdummy�����ӵ�λ��
	vector<vec3> uArray;					//��¼kʱ���������ӵ��ٶ�
	vector<vec3> uArray1;					//��¼kʱ�̳�ȥdummy�����ӵ��ٶ�
	vector<vec3> tempUArray;				//��¼����ѹ�����k+1ʱ���������ӵ�u*
	vector<vec3> tempUArray1;				//��¼����ѹ�����k+1ʱ�̳�ȥdummy�����ӵ�u*
	vector<vec3> tempPosArray;				//��¼����ѹ�����k+1ʱ���������ӵ�r*
	vector<vec3> tempPosArray1;				//��¼����ѹ�����k+1ʱ�̳�ȥdummy�����ӵ�r*

	//1.2 ������ʽ��P
	vector<bool> surfaceJudgeArray;			//��¼���ӵı����ж�
	vector<float> n0Array;					//��¼���ӵĳ�ʼ�ܶ�

	for (int i = 0; i < particles.size(); i++)
	{
		posArray.push_back(particles[i].position);
		uArray.push_back(particles[i].speed);
		if (i < (particleNumber + particleWallNum))
		{
			posArray1.push_back(particles[i].position);
			uArray1.push_back(particles[i].speed);
		}
	}
	//�������� ��ȡ��ʱ���ٶȺ�λ�ã�ֻ����ͨ������Ҫ���㣬�Ҽ����ʱ����Ҫdummy���������,��wall���ӵ��ٶȺ�λ�ö����䣩
	for (int i = 0; i < (particleNumber + particleWallNum); i++)
	{
		if (i < particleNumber)
		{
			//������ʱ�ٶ�u* ��������ʱλ��r*
			vec3 tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(uArray1, posArray1, i, particles[i].n0), particles[i].position);
			tempUArray1.push_back(tempU);

			vec3 tempPos = particles[i].position + tempU * mpsTool->GetDeltaT();
			tempPosArray1.push_back(tempPos);
		}
		else
		{
			tempUArray1.push_back(vec3(0));						//wall�����ٶ�һֱΪ0
			tempPosArray1.push_back(particles[i].position);		//wall����λ�ñ��ֲ���
		}
	}
	//Ϊȫ��������ʱֵ���ϸ�ֵ
	tempUArray = tempUArray1;
	tempPosArray = tempPosArray1;
	for (int i = (particleNumber + particleWallNum); i < particleDummyNum; i++)
	{
		tempUArray.push_back(vec3(0));
		tempPosArray.push_back(particles[i].position);
	}


	//����ʱ�ٶȺ���ʱλ�ü���ÿ�����Ӷ�Ӧ���Ҷ���
	for (int i = 0; i < (particleNumber + particleWallNum); i++)
	{
		//��ͨ���Ӻ�wall���ӵ��Ҷ��������Ҫ��ͬ�����Ӽ���
		if (i < particleNumber)
		{
			//��ͨ���Ӳ���Ҫdummy
			float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray1, i));
			Right.push_back(resRight);
		}
		else
		{
			//wall ������ȫ�����Ӽ���
			float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
			Right.push_back(resRight);
		}
		//��ʼ����ʱ���þɵ��Ҷ����
		float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
		//float resRight = mpsTool->ImplicitLaplacianRight(particles[i].n0, resDivergence, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
		Right.push_back(resRight);

		//��Ϊ���������Ҷ���ļ��㲻Ӱ�죬����ͬһ��ѭ�������Ч�ʣ���ʼ����ʱ����Ҫ�ٽ���һ�α�����
		surfaceJudgeArray.push_back(particles[i].isSurface);
		n0Array.push_back(particles[i].n0);
	}

	//1.2.2 ��һ�����ɷ���(��ʱ�ļ����ų�dummy)
	vector<double> resP = mpsTool->ImplicitCalculateP(posArray1, n0Array, surfaceJudgeArray, Right);

	for (int i = 0; i < (particleNumber + particleWallNum); i++)
	{
		particles[i].pressure = resP[i];
	}
}

void MPSWaterParticleGroup::InitMPSTool()
{
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	mpsTool->SetRe(l0);
	mpsTool->SetViscosity(viscosity);
	mpsTool->SetDeltaT(0.01);
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
	marchingCube = new MarchingCube();
	marchingCube->GetMeshData(posArray, verticesInfo, provideNormal, provideTex);
	//MarchingCube::GetInstance()->GetMeshData(posArray, verticesInfo, provideNormal, provideTex);

	shaderData->InitVertexBuffer(verticesInfo, provideNormal, provideTex);
}

void MPSWaterParticleGroup::SetInitialN0()
{
	MPSToolFun* tool = MPSToolFun::GetMPSTool();
	vector<vec3> posArray;								//�������ӵ�λ�ü��ϣ���Ҫ��wall������Ҫ
	vector<vec3> posArray1;								//ȥ��dummy���ӵļ��ϣ�������ͨ����
	for (int i = 0; i < particles.size(); i++)
	{
		posArray.push_back(particles[i].position);
		if (i < (particleNumber + particleWallNum))
			posArray1.push_back(particles[i].position);
	}

	//ֻ��Ҫ������ͨ���Ӻ�ǽ���ӵ��ܶ�
	for (int i = 0; i < particles.size(); i++)
	{
		if (i < particleNumber)			//��ͨ����
			particles[i].n0 = tool->DensityN(posArray1, i);
		else if (i < (particleNumber + particleWallNum))		//ǽ����
			particles[i].n0 = tool->DensityN(posArray, i);
	}

}

void MPSWaterParticleGroup::UpdateAdjoin(float range)
{
	int num1 = particleNumber + particleWallNum;
	int num2 = particles.size();
	//�ڽ���������ͨ���Ӳ���Ҫ����dummy��wall���Ӷ�Ҫ����
	for (int i = 0; i < num1; i++)
	{
		particles[i].adjoinParticleIndex.clear();
		if (i < particleNumber)
		{
			//��ͨ���ӵĴ���
			for (int j = 0; j < num1; j++)
			{
				if (j != i)
				{
					vec3 pos1 = particles[j].position;
					if (distance(pos1, particles[i].position) <= range)
						particles[i].adjoinParticleIndex.push_back(j);
				}
			}
		}
		else
		{
			//wall���ӵĴ���
			for (int j = 0; j < num2; j++)
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

}

void MPSWaterParticleGroup::Update(float dt)
{
	//ÿһ֡�ļ�����Ҫ�ȸ���һ���ڽӹ�ϵ
	UpdateAdjoin(range);

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
	vector<vec3> tempUArray;				//��¼����ѹ�����k+1ʱ���������ӵ�u*
	vector<vec3> tempPosArray;				//��¼����ѹ�����k+1ʱ���������ӵ�r*

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
		vec3 tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(uArray, posArray, i, particles[i].n0), particles[i].position);
		tempUArray.push_back(tempU);

		vec3 tempPos = particles[i].position + tempU * mpsTool->GetDeltaT();;
		tempPosArray.push_back(tempPos);

		//����ʱ�̵�λ�ø���һ�����ӵ��ܶ�
		particles[i].tho = mpsTool->DensityN(posArray, i);
	}
	//����ʱ�ٶȺ���ʱλ�ü���ÿ�����Ӷ�Ӧ���Ҷ���
	for (int i = 0; i < particles.size(); i++)
	{
		float resDivergence = mpsTool->ExplicitDivergence(tempUArray, tempPosArray, i, particles[i].n0);
		float resRight = mpsTool->ImplicitLaplacianRight(particles[i].n0, resDivergence, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
		Right.push_back(resRight);

		//��Ϊ���������Ҷ���ļ��㲻Ӱ�죬����ͬһ��ѭ�������Ч��
		particles[i].SurfaceAdjudge(a, g, l0);
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
		particles[i].speed = mpsTool->CalculateU(resLU, v1, particles[i].speed, particles[i].tho);
	}

	//����λ�ú�ѹ��
	for (int i = 0; i < particles.size(); i++)
	{
		particles[i].pressure = resP[i];
		particles[i].position += particles[i].speed * mpsTool->GetDeltaT();
	}

	//������õ�λ�õ���н�ģ
	//Modeling();
}

void MPSWaterParticleGroup::Draw()
{

}

void MPSWaterParticleGroup::InitBufferData()
{
}
