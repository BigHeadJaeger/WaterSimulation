#include "ParticleGroup.h"

void MPSWaterParticleGroup::InitParticles()
{
	//��ʼ��MPS����
	InitMPSTool();
	ParticlesArrange();
	SetInitialN0();

	//�����ڽӵ��ϵ
	//UpdateAdjoin(range);
	//����ÿһ�����ӵĳ�ʼ�ܶ�
	//SetInitialN0();


	////�˴������������ӵĳ�ʼѹ��
	//MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	//vector<double> Right;					//�洢ÿһ�����ӵ��Ҷ���

	//vector<vec3> posArray;					//��¼kʱ���������ӵ�λ��
	//vector<vec3> posArray1;					//��¼kʱ�̳�ȥdummy�����ӵ�λ��
	//vector<vec3> uArray;					//��¼kʱ���������ӵ��ٶ�
	//vector<vec3> uArray1;					//��¼kʱ�̳�ȥdummy�����ӵ��ٶ�
	//vector<vec3> tempUArray;				//��¼����ѹ�����k+1ʱ���������ӵ�u*
	//vector<vec3> tempUArray1;				//��¼����ѹ�����k+1ʱ�̳�ȥdummy�����ӵ�u*
	//vector<vec3> tempPosArray;				//��¼����ѹ�����k+1ʱ���������ӵ�r*
	//vector<vec3> tempPosArray1;				//��¼����ѹ�����k+1ʱ�̳�ȥdummy�����ӵ�r*

	////1.2 ������ʽ��P
	//vector<bool> surfaceJudgeArray;			//��¼���ӵı����ж�
	//vector<float> n0Array;					//��¼���ӵĳ�ʼ�ܶ�

	//for (int i = 0; i < particles.size(); i++)
	//{
	//	posArray.push_back(particles[i].position);
	//	uArray.push_back(particles[i].speed);
	//	if (i < (particleNumber + particleWallNum))
	//	{
	//		posArray1.push_back(particles[i].position);
	//		uArray1.push_back(particles[i].speed);
	//	}
	//}
	////�������� ��ȡ��ʱ���ٶȺ�λ�ã�ֻ����ͨ������Ҫ���㣬�Ҽ����ʱ����Ҫdummy���������,��wall���ӵ��ٶȺ�λ�ö����䣩
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	if (i < particleNumber)
	//	{
	//		//������ʱ�ٶ�u* ��������ʱλ��r*
	//		vec3 tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(uArray1, posArray1, i, particles[i].n0), particles[i].speed);
	//		tempUArray1.push_back(tempU);

	//		vec3 tempPos = particles[i].position + tempU * mpsTool->GetDeltaT();
	//		tempPosArray1.push_back(tempPos);
	//	}
	//	else
	//	{
	//		tempUArray1.push_back(vec3(0));						//wall�����ٶ�һֱΪ0
	//		tempPosArray1.push_back(particles[i].position);		//wall����λ�ñ��ֲ���
	//	}
	//}
	////Ϊȫ��������ʱֵ���ϸ�ֵ
	//tempUArray = tempUArray1;
	//tempPosArray = tempPosArray1;
	//for (int i = (particleNumber + particleWallNum); i < particleTotalNum; i++)
	//{
	//	tempUArray.push_back(vec3(0));
	//	tempPosArray.push_back(particles[i].position);
	//}


	////����ʱ�ٶȺ���ʱλ�ü���ÿ�����Ӷ�Ӧ���Ҷ���
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	//��ͨ���Ӻ�wall���ӵ��Ҷ��������Ҫ��ͬ�����Ӽ���
	//	//��ʼ����ʱ���þɵ��Ҷ����
	//	if (i < particleNumber)
	//	{
	//		//��ͨ���Ӳ���Ҫdummy
	//		float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray1, i));
	//		Right.push_back(resRight);
	//	}
	//	else
	//	{
	//		//wall ������ȫ�����Ӽ���
	//		float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
	//		Right.push_back(resRight);
	//	}
	//	
	//	//��Ϊ���������Ҷ���ļ��㲻Ӱ�죬����ͬһ��ѭ�������Ч�ʣ���ʼ����ʱ����Ҫ�ٽ���һ�α�����
	//	surfaceJudgeArray.push_back(particles[i].isSurface);
	//	n0Array.push_back(particles[i].n0);
	//}

	////1.2.2 ��һ�����ɷ���(��ʱ�ļ����ų�dummy)
	//vector<double> resP = mpsTool->ImplicitCalculateP(tempPosArray1, n0Array, surfaceJudgeArray, Right);
	////vector<double> resP1 = mpsTool->ImplicitCalculateP(tempPosArray1, n0Array, surfaceJudgeArray, Right);
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	if (particles[i].isSurface)
	//		particles[i].pressure = 0;
	//	else
	//		particles[i].pressure = resP[i];
	//}
}

void MPSWaterParticleGroup::ParticlesArrange()
{
	vec3 initPos = transformation.position;
	vec3 offset = vec3(l0);									//���ӿռ�ֲ�����
	int width, height, depth;								//���ӵĳ�����Ų�
	width = height = depth = 4;

	height = 5;

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
			R.particleType = SURFACE;
		else
			R.particleType = FLUID;
		R.index = i;
		particles.push_back(R);
	}
	particleNumber = n;

	//����Ҫ��ʼ������ǽ�Լ�dummy wall	(���Կ����ٴ���1+2������޸ǵ�cube)   ps:Ҫȷ��������Χס���е�����
	//wall
	int containerWidth = width + 2;
	int containerHeight = height + 3;
	int containerDepth = depth + 2;
	initPos -= offset;
	positions.clear();
	n = UncoverHollowCubeDistribute(positions, initPos, offset, containerWidth, containerHeight, containerDepth);
	for (int i = 0; i < n; i++)
	{
		MPSWaterParticle R;
		R.position = positions[i];
		R.particleType = WALL;
		R.index = i + particleNumber;				//����֮ǰ�Ѿ��е�����ƫ����
		particles.push_back(R);
	}
	particleWallNum = n;

	containerWidth += 2;
	containerHeight += 2;
	containerDepth += 2;
	initPos -= offset;
	positions.clear();
	n = UncoverHollowCubeDistribute(positions, initPos, offset, containerWidth, containerHeight, containerDepth);
	for (int i = 0; i < n; i++)
	{
		MPSWaterParticle R;
		R.position = positions[i];
		R.particleType = WALL;
		R.index = i + particleNumber;				//����֮ǰ�Ѿ��е�����ƫ����
		particles.push_back(R);
	}
	particleWallNum += n;

	//dummy
	int dummyWidth = containerWidth + 2;
	int dummyHeight = containerHeight + 2;
	int dummyDepth = containerDepth + 2;
	positions.clear();
	initPos -= offset;

	n = UncoverHollowCubeDistribute(positions, initPos, offset, dummyWidth, dummyHeight, dummyDepth);
	//�ڶ���
	dummyHeight += 1;
	dummyWidth += 2;
	dummyDepth += 2;
	initPos -= offset;

	n += UncoverHollowCubeDistribute(positions, initPos, offset, dummyWidth, dummyHeight, dummyDepth);

	//dummyHeight += 1;
	//dummyWidth += 2;
	//dummyDepth += 2;
	//initPos -= offset;

	//n += UncoverHollowCubeDistribute(positions, initPos, offset, dummyWidth, dummyHeight, dummyDepth);
	for (int i = 0; i < n; i++)
	{
		MPSWaterParticle R;
		R.position = positions[i];
		R.particleType = DUMMY;
		R.index = i + particleNumber + particleWallNum;				//����֮ǰ�Ѿ��е�����ƫ����
		particles.push_back(R);
	}

	particleDummyNum = n;

	particleTotalNum = particleNumber + particleDummyNum + particleWallNum;
}

void MPSWaterParticleGroup::InitMPSTool()
{
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
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

// �����Ӧ��3��n0ֵ�����㷽��Ϊ�����һ������ȫ��Χ���м�����ӣ�����Ӱ�췶Χ�ڳ������ӣ�����������ӵ��ܶ���Ϊ��ʼ�ܶ�
void MPSWaterParticleGroup::SetInitialN0()
{
	MPSToolFun* tool = MPSToolFun::GetMPSTool();

	vec3 ri(0, 0, 0), rj(0, 0, 0);
	vec3 r_min(-4);
	vec3 r_max(5);
	double dis;

	if (DIMENSION == 2) {
		r_min.z = 0;
		r_max.z = 1;
	}

	//r_min.z = 0;
	//r_max.z = 1;

	for (int x = r_min.x; x < r_max.x; ++x)
	{
		for (int y = r_min.y; y < r_max.y; ++y)
		{
			for (int z = r_min.z; z < r_max.z; ++z)
			{
				if (((x == 0) && (y == 0)) && (z == 0))continue;
				rj = (float)l0 * vec3(x, y, z);
				dis = distance(rj, ri);
				
				n0ForNumberDensity += tool->WeightFun(dis, RADIUS_FOR_NUMBER_DENSITY);
				n0ForGradient += tool->WeightFun(dis, RADIUS_FOR_GRADIENT);
				//n0ForLambda += dis * dis * tool->WeightFun(dis, RADIUS_FOR_LAPLACIAN);
				n0ForLambda += tool->WeightFun(dis, RADIUS_FOR_LAPLACIAN);

				lambda0 += dis * dis * tool->WeightFun(dis, RADIUS_FOR_LAPLACIAN);
			}
		}
	}

	lambda0 /= n0ForLambda;

	tool->testLambda0 = lambda0;



	//vector<vec3> posArray;								//�������ӵ�λ�ü��ϣ���Ҫ��wall������Ҫ
	//vector<vec3> posArray1;								//ȥ��dummy���ӵļ��ϣ�������ͨ����
	//for (int i = 0; i < particles.size(); i++)
	//{
	//	posArray.push_back(particles[i].position);
	//	if (i < (particleNumber + particleWallNum))
	//		posArray1.push_back(particles[i].position);
	//}

	////ֻ��Ҫ������ͨ���Ӻ�ǽ���ӵ��ܶ�
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	if (i < particleNumber)			//��ͨ����
	//		particles[i].n0 = tool->DensityN(posArray1, i);
	//	else							//ǽ����
	//		particles[i].n0 = tool->DensityN(posArray, i);
	//}

}

//void MPSWaterParticleGroup::UpdateAdjoin(float range)
//{
//	int num1 = particleNumber + particleWallNum;
//	int num2 = particles.size();
//	//�ڽ���������ͨ���Ӳ���Ҫ����dummy��wall���Ӷ�Ҫ����
//	for (int i = 0; i < num1; i++)
//	{
//		particles[i].adjoinParticleIndex.clear();
//		if (i < particleNumber)
//		{
//			//��ͨ���ӵĴ���
//			for (int j = 0; j < num1; j++)
//			{
//				if (j != i)
//				{
//					vec3 pos1 = particles[j].position;
//					if (distance(pos1, particles[i].position) <= range)
//						particles[i].adjoinParticleIndex.push_back(j);
//				}
//			}
//		}
//		else
//		{
//			//wall���ӵĴ���
//			for (int j = 0; j < num2; j++)
//			{
//				if (j != i)
//				{
//					vec3 pos1 = particles[j].position;
//					if (distance(pos1, particles[i].position) <= range)
//						particles[i].adjoinParticleIndex.push_back(j);
//				}
//			}
//		}
//	}
//
//}

void MPSWaterParticleGroup::Update(float dt)
{
	shaderData->UpdateMatrix(transformation);
	//return;
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	//������ʾ���裨����ѹ����ĵ�һ���ٶ�����У����

	vector<vec3> initUArrayWallAndFluid;			// ʱ�䲽n���ٶȼ���
	vector<vec3> initPArrayWallAndFluid;			// ʱ�䲽n��λ�ü���
	vector<vec3> middleUArray;						//
	vector<vec3> middlePArray;						//

	vector<vec3> middlePAllArray;					// �������ӵ��м�ֵ���꣨��������Wall���м��ܶȣ�

	vector<double> Right;					//�洢ÿһ�����ӵ��Ҷ���
	vector<bool> ignoreArray;				// ѹ�����ɷ�������Ҫ���Ե��surface��

	for (int i = 0; i < particleNumber + particleWallNum; i++)
	{
		initUArrayWallAndFluid.push_back(particles[i].speed);
		initPArrayWallAndFluid.push_back(particles[i].position);
	}

	// ����ÿ���������ӵ��м��ٶ�
	for (int i = 0; i < particleNumber + particleWallNum; i++)
	{
		vec3 tempU = vec3(particles[i].speed);
		vec3 tempPos = vec3(particles[i].position);
		if (i< particleNumber)
		{
			tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(initUArrayWallAndFluid, initPArrayWallAndFluid, i, n0ForLambda), particles[i].speed);
			//tempPos = particles[i].position + (float)0.5 * (tempU + particles[i].speed) * (float)DELTA_TIME;
			tempPos = particles[i].position + tempU * (float)DELTA_TIME;
		}


		middleUArray.push_back(tempU);
		middlePArray.push_back(tempPos);

		middlePAllArray.push_back(tempPos);
	}

	// ��dummy����
	for (int i = particleNumber + particleWallNum; i < particleTotalNum; i++)
		middlePAllArray.push_back(particles[i].position);

	// Fluid���м��������ܶ�
	for (int i = 0; i < particleNumber; i++)
	{
		float tempN = mpsTool->DensityN(middlePArray, i);
		particles[i].middleN = tempN;
		//cout << n0ForNumberDensity << " " << tempN << endl;
		// ������
		if (tempN < n0ForNumberDensity * THRESHOLD_RATIO_OF_NUMBER_DENSITY)
		{
			particles[i].particleType = SURFACE;
			ignoreArray.push_back(true);
		}
		else
		{
			particles[i].particleType == FLUID;
			ignoreArray.push_back(false);
		}

			


	}
	// Wall���м��������ܶ�
	for (int i = particleNumber; i < particleNumber + particleWallNum; i++)
	{
		float tempN = mpsTool->DensityN(middlePAllArray, i);
		particles[i].middleN = tempN;
		ignoreArray.push_back(false);
		//if (tempN < n0ForNumberDensity * THRESHOLD_RATIO_OF_NUMBER_DENSITY)
		//{
		//	particles[i].particleType = SURFACE;
		//	ignoreArray.push_back(true);
		//}
		//else
		//{
		//	particles[i].particleType == FLUID;
		//	ignoreArray.push_back(false);
		//}
	}


	//for (int i = 0; i < particleNumber + particleWallNum; i++)
	//{
	//	particles[i].speed = middleUArray[i];
	//	particles[i].position = middlePArray[i];
	//}

	//InitBufferData();
	//return;


	// Wall��Fluid���ӵ��Ҷ������Ҫ�ų���surface�ļ��㣩
	for (int i = 0; i < particleNumber + particleWallNum; i++)
	{
		float resRight = 0;
		//if (!ignoreArray[i])
		//{

		//}
		resRight = mpsTool->OldImplicitLaplacianRight(FLUID_DENSITY, n0ForNumberDensity, particles[i].middleN);
		Right.push_back(resRight);
	}



	// ��һ�����ɷ���
	vector<double> resP = mpsTool->ImplicitCalculateP(middlePArray, n0ForLambda, ignoreArray, Right);

	for (int i = 0; i < resP.size(); i++)
	{
		if (ignoreArray[i])
			resP[i] = 0;
	}


	for (int i = 0; i < particleNumber + particleWallNum; i++)
	{
		vec3 secondTempU = middleUArray[i];
		vec3 secondTempP = middlePArray[i];
		if (i < particleNumber)
		{
			vec3 resGradient = mpsTool->OldGradient(middlePArray, resP, i, n0ForGradient);

			secondTempU = (float)(-DELTA_TIME / FLUID_DENSITY) * resGradient + middleUArray[i];
			//secondTempP = middlePArray[i] + (float)0.5 * (middleUArray[i] + secondTempU) * (float)DELTA_TIME;
			//secondTempP = middlePArray[i] + secondTempU * (float)DELTA_TIME;
			secondTempP = middlePArray[i] + (secondTempU - middleUArray[i]) * (float)DELTA_TIME;
		}

		particles[i].speed = secondTempU;
		particles[i].position = secondTempP;

	}

	printf("sds");






	//MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	//ÿһ֡����Ҫ�㷨���̱���ÿһ������

	//1.����ÿ�����ӵ�p
	//1.1����ÿ�����ӵ��Ҷ���
	// u*��ɢ��
	// u��������˹���
	// u*�ļ���
	// n*�ļ���
	//MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	//vector<double> Right;					//�洢ÿһ�����ӵ��Ҷ���

	//vector<vec3> posArray;					//��¼kʱ���������ӵ�λ��
	//vector<vec3> posArray1;					//��¼kʱ�̳�ȥdummy�����ӵ�λ��
	//vector<vec3> uArray;					//��¼kʱ���������ӵ��ٶ�
	//vector<vec3> uArray1;					//��¼kʱ�̳�ȥdummy�����ӵ��ٶ�
	//vector<vec3> tempUArray;				//��¼����ѹ�����k+1ʱ���������ӵ�u*
	//vector<vec3> tempUArray1;				//��¼����ѹ�����k+1ʱ�̳�ȥdummy�����ӵ�u*
	//vector<vec3> tempPosArray;				//��¼����ѹ�����k+1ʱ���������ӵ�r*
	//vector<vec3> tempPosArray1;				//��¼����ѹ�����k+1ʱ�̳�ȥdummy�����ӵ�r*

	////1.2 ������ʽ��P
	//vector<bool> surfaceJudgeArray;			//��¼���ӵı����ж�
	//vector<float> n0Array;					//��¼���ӵĳ�ʼ�ܶ�

	//for (int i = 0; i < particles.size(); i++)
	//{
	//	posArray.push_back(particles[i].position);
	//	uArray.push_back(particles[i].speed);
	//	if (i < (particleNumber + particleWallNum))
	//	{
	//		posArray1.push_back(particles[i].position);
	//		uArray1.push_back(particles[i].speed);
	//	}
	//}
	////�������� ��ȡ��ʱ���ٶȺ�λ�ã�ֻ����ͨ������Ҫ���㣬�Ҽ����ʱ����Ҫdummy���������,��wall���ӵ��ٶȺ�λ�ö����䣩
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	//������ʱֵ��ʱ��˳�����һ����ʱ�̵������ܶ�
	//	if (i < particleNumber)
	//	{
	//		//������ʱ�ٶ�u* ��������ʱλ��r*
	//		vec3 tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(uArray1, posArray1, i, particles[i].n0), particles[i].speed);
	//		tempUArray1.push_back(tempU);

	//		vec3 tempPos = particles[i].position + tempU * mpsTool->GetDeltaT();
	//		tempPosArray1.push_back(tempPos);

	//		particles[i].tho = mpsTool->DensityN(posArray1, i);
	//	}
	//	else
	//	{
	//		tempUArray1.push_back(vec3(0));						//wall�����ٶ�һֱΪ0
	//		tempPosArray1.push_back(particles[i].position);		//wall����λ�ñ��ֲ���
	//		particles[i].tho = mpsTool->DensityN(posArray, i);
	//	}
	//}
	////Ϊȫ��������ʱֵ���ϸ�ֵ
	//tempUArray = tempUArray1;
	//tempPosArray = tempPosArray1;
	//for (int i = (particleNumber + particleWallNum); i < particleTotalNum; i++)
	//{
	//	tempUArray.push_back(vec3(0));
	//	tempPosArray.push_back(particles[i].position);
	//}

	////vector<float> tempNArray;
	////����ʱ�ٶȺ���ʱλ�ü���ÿ�����Ӷ�Ӧ���Ҷ���
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	float tempN = 0;
	//	//��ͨ���Ӻ�wall���ӵ��Ҷ��������Ҫ��ͬ�����Ӽ���
	//	if (i < particleNumber)
	//	{
	//		//��ͨ���Ӳ���Ҫdummy
	//		float resDivergence = mpsTool->ExplicitDivergence(tempUArray1, tempPosArray1, i, particles[i].n0);
	//		tempN = mpsTool->DensityN(tempPosArray1, i);
	//		//tempNArray.push_back(tempN);
	//		float resRight = mpsTool->ImplicitLaplacianRight(particles[i].n0, resDivergence, particles[i].n0, tempN);
	//		Right.push_back(resRight);
	//	}
	//	else
	//	{
	//		//wall ������ȫ�����Ӽ���
	//		float resDivergence = mpsTool->ExplicitDivergence(tempUArray, tempPosArray, i, particles[i].n0);
	//		tempN = mpsTool->DensityN(tempPosArray, i);
	//		//tempNArray.push_back(tempN);
	//		float resRight = mpsTool->ImplicitLaplacianRight(particles[i].n0, resDivergence, particles[i].n0, tempN);
	//		Right.push_back(resRight);
	//	}

	//	//��Ϊ���������Ҷ���ļ��㲻Ӱ�죬����ͬһ��ѭ�������Ч��
	//	//����ʱ��wallҲ��⣬��Ϊ��dummy���Բ��ᱻ����Ϊ�Ǳ���
	//	
	//	//particles[i].OldSurfaceAdjudge(0.97, tempN);
	//	if (i < particleNumber)
	//	{
	//		//particles[i].SurfaceAdjudge(a, g, l0, particles[i].n0);
	//		if (particles[i].particleType == SURFACE)
	//			surfaceJudgeArray.push_back(true);
	//	}
	//	else
	//		surfaceJudgeArray.push_back(false);
	//	n0Array.push_back(particles[i].n0);
	//}

	////1.2.2 ��һ�����ɷ���(��ʱ�ļ����ų�dummy)
	//vector<double> resP = mpsTool->ImplicitCalculateP(posArray1, n0Array, surfaceJudgeArray, Right);

	////����ÿһ�������µ��ٶ�U��ֻ��Ҫ������ͨ���ӾͿ��ԣ�����wall������Ҫ�������ʽ��
	//for (int i = 0; i < particleNumber; i++)
	//{
	//	mat3 C = mpsTool->GetMaterixC(posArray1, i, particles[i].n0);
	//	vec3 v1 = mpsTool->ExplicitGradient(C, resP, posArray1, particles[i].n0, i);			//����ÿ�����ӵ���ʽ�ݶ�
	//	vec3 resLU = mpsTool->ExplicitLaplacian(uArray1, posArray1, i, particles[i].n0);
	//	particles[i].speed = mpsTool->CalculateU(resLU, v1, particles[i].speed, particles[i].tho);
	//}

	////����λ�ú�ѹ����ֻ����ͨ���ӵ�λ�ú�ѹ���б仯��wall����ֻ��ѹ����䣩
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	if (particles[i].particleType == SURFACE)
	//		particles[i].pressure = 0;
	//	else
	//		particles[i].pressure = resP[i];
	//	particles[i].position += particles[i].speed * mpsTool->GetDeltaT();
	//}

	//������õ�λ�õ���н�ģ
	//Modeling();
	InitBufferData();
}

void MPSWaterParticleGroup::Draw()
{
	// GL_POINTS
	shaderData->SetDrawType(GL_POINTS);
	renderer->Render(shaderData);
}

void MPSWaterParticleGroup::InitBufferData()
{
	vector<float> data;

	for (int i = 0; i < particles.size(); i++)
	{
		vec3 color = vec3(0);
		switch (particles[i].particleType)
		{
		case FLUID:
			color = vec3(0, 0, 255);
			break;
		case SURFACE:
			color = vec3(255, 0, 0);
			break;
		case DUMMY:
			color = vec3(0);
			//color = vec3(0, 255, 255);
			break;
		case WALL:
			color = vec3(0, 255, 0);
			break;
		default:
			break;
		}

		data.push_back(particles[i].position.x);
		data.push_back(particles[i].position.y);
		data.push_back(particles[i].position.z);
		data.push_back(color.x);
		data.push_back(color.y);
		data.push_back(color.z);
	}

	shaderData->drawUnitNumber = data.size() / 6;
	dynamic_cast<VCShaderData*>(shaderData)->InitVertexBuffer(data, GL_DYNAMIC_DRAW);
}

void MPSWaterParticleGroup::UpdateBufferData()
{
}
