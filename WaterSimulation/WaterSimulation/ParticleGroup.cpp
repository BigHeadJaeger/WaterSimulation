#include "ParticleGroup.h"

void MPSWaterParticleGroup::InitParticles()
{
	//初始化MPS工具
	InitMPSTool();
	vec3 initPos = transformation.position;
	vec3 offset = vec3(l0);									//粒子空间分布布局
	int width, height, depth;								//粒子的长宽高排布
	width = height = depth = 10;

	float highest = 0;										//找到最高点位置
	if (offset.y > 0)
		highest = initPos.y + (height - 1) * offset.y;
	else
		highest = initPos.y;
	vector<vec3> positions;
	int n = 0;
	n = CubeDistribute(positions, initPos, offset, width, height, depth);
	for (int i = 0; i < n; i++)
	{
		//初始化当前粒子并push进去
		MPSWaterParticle R;
		//此处指定粒子属性 位置
		R.position = positions[i];
		//如果这个粒子在最上面则标记为表面粒子
		if (positions[i].y == highest)
			R.isSurface = true;
		R.index = i;
		particles.push_back(R);
	}
	particleNumber = n;

	//还需要初始化粒子墙以及dummy wall	(可以看成再创建1+2层空心无盖的cube)   ps:要确保容器包围住所有的粒子
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
		R.index = i + particleNumber;				//加上之前已经有的索引偏移量
	}
	particleWallNum = n;
	//dummy
	int dummyWidth = containerWidth + 2;
	int dummyHeight = containerHeight;
	int dummyDepth = containerDepth + 2;
	positions.clear();

	initPos -= offset;
	n = UncoverHollowCubeDistribute(positions, initPos, offset, containerWidth, containerHeight, containerDepth);
	//第二层
	dummyWidth += 2;
	dummyDepth += 2;
	initPos -= offset;
	n += UncoverHollowCubeDistribute(positions, initPos, offset, containerWidth, containerHeight, containerDepth);
	for (int i = 0; i < n; i++)
	{
		MPSWaterParticle R;
		R.position = positions[i];
		R.isDummy = true;
		R.index = i + particleNumber + particleWallNum;				//加上之前已经有的索引偏移量
	}
	particleDummyNum = n;

	particleTotalNum = particleNumber + particleDummyNum + particleWallNum;

	//更新邻接点关系
	UpdateAdjoin(range);
	//设置每一个粒子的初始密度
	SetInitialN0();


	//此处计算所有粒子的初始压力
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	vector<double> Right;					//存储每一个粒子的右端项

	vector<vec3> posArray;					//记录k时刻所有粒子的位置
	vector<vec3> posArray1;					//记录k时刻除去dummy的粒子的位置
	vector<vec3> uArray;					//记录k时刻所有粒子的速度
	vector<vec3> uArray1;					//记录k时刻除去dummy的粒子的速度
	vector<vec3> tempUArray;				//记录忽略压力项后k+1时刻所有粒子的u*
	vector<vec3> tempUArray1;				//记录忽略压力项后k+1时刻除去dummy的粒子的u*
	vector<vec3> tempPosArray;				//记录忽略压力项后k+1时刻所有粒子的r*
	vector<vec3> tempPosArray1;				//记录忽略压力项后k+1时刻除去dummy的粒子的r*

	//1.2 计算隐式的P
	vector<bool> surfaceJudgeArray;			//记录粒子的表面判断
	vector<float> n0Array;					//记录粒子的初始密度

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
	//遍历粒子 获取临时的速度和位置（只有普通粒子需要计算，且计算的时候不需要dummy的相关数据,而wall粒子的速度和位置都不变）
	for (int i = 0; i < (particleNumber + particleWallNum); i++)
	{
		if (i < particleNumber)
		{
			//计算临时速度u* 并更新临时位置r*
			vec3 tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(uArray1, posArray1, i, particles[i].n0), particles[i].position);
			tempUArray1.push_back(tempU);

			vec3 tempPos = particles[i].position + tempU * mpsTool->GetDeltaT();
			tempPosArray1.push_back(tempPos);
		}
		else
		{
			tempUArray1.push_back(vec3(0));						//wall粒子速度一直为0
			tempPosArray1.push_back(particles[i].position);		//wall粒子位置保持不变
		}
	}
	//为全部粒子临时值集合赋值
	tempUArray = tempUArray1;
	tempPosArray = tempPosArray1;
	for (int i = (particleNumber + particleWallNum); i < particleDummyNum; i++)
	{
		tempUArray.push_back(vec3(0));
		tempPosArray.push_back(particles[i].position);
	}


	//用临时速度和临时位置计算每个粒子对应的右端项
	for (int i = 0; i < (particleNumber + particleWallNum); i++)
	{
		//普通粒子和wall粒子的右端项计算需要不同的粒子集合
		if (i < particleNumber)
		{
			//普通粒子不需要dummy
			float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray1, i));
			Right.push_back(resRight);
		}
		else
		{
			//wall 粒子用全部粒子集合
			float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
			Right.push_back(resRight);
		}
		//初始化的时候用旧的右端项方法
		float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
		//float resRight = mpsTool->ImplicitLaplacianRight(particles[i].n0, resDivergence, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
		Right.push_back(resRight);

		//因为表面检测与右端项的计算不影响，放在同一个循环中提高效率，初始化的时候不需要再进行一次表面检测
		surfaceJudgeArray.push_back(particles[i].isSurface);
		n0Array.push_back(particles[i].n0);
	}

	//1.2.2 解一个泊松方程(此时的计算排除dummy)
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
	vector<vec3> posArray;								//所有粒子的位置集合，主要是wall粒子需要
	vector<vec3> posArray1;								//去除dummy粒子的集合，用于普通粒子
	for (int i = 0; i < particles.size(); i++)
	{
		posArray.push_back(particles[i].position);
		if (i < (particleNumber + particleWallNum))
			posArray1.push_back(particles[i].position);
	}

	//只需要计算普通粒子和墙粒子的密度
	for (int i = 0; i < particles.size(); i++)
	{
		if (i < particleNumber)			//普通粒子
			particles[i].n0 = tool->DensityN(posArray1, i);
		else if (i < (particleNumber + particleWallNum))		//墙粒子
			particles[i].n0 = tool->DensityN(posArray, i);
	}

}

void MPSWaterParticleGroup::UpdateAdjoin(float range)
{
	int num1 = particleNumber + particleWallNum;
	int num2 = particles.size();
	//邻接粒子中普通粒子不需要考虑dummy，wall粒子都要考虑
	for (int i = 0; i < num1; i++)
	{
		particles[i].adjoinParticleIndex.clear();
		if (i < particleNumber)
		{
			//普通粒子的处理
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
			//wall粒子的处理
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
	//每一帧的计算需要先更新一下邻接关系
	UpdateAdjoin(range);

	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	//每一帧的主要算法流程遍历每一个粒子

	//1.计算每个粒子的p
	//1.1计算每个粒子的右端项
	// u*的散度
	// u的拉普拉斯结果
	// u*的计算
	// n*的计算
	vector<double> Right;					//存储每一个粒子的右端项
	
	vector<vec3> posArray;					//记录k时刻所有粒子的位置
	vector<vec3> uArray;					//记录k时刻所有粒子的速度
	vector<vec3> tempUArray;				//记录忽略压力项后k+1时刻所有粒子的u*
	vector<vec3> tempPosArray;				//记录忽略压力项后k+1时刻所有粒子的r*

	//1.2 计算隐式的P
	vector<bool> surfaceJudgeArray;			//记录粒子的表面判断
	vector<float> n0Array;					//记录粒子的初始密度

	for (int i = 0; i < particles.size(); i++)
	{
		posArray.push_back(particles[i].position);
		uArray.push_back(particles[i].speed);
	}
	//遍历粒子 获取临时的速度和位置
	for (int i = 0; i < particles.size(); i++)
	{
		//计算临时速度u* 并更新临时位置r*
		vec3 tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(uArray, posArray, i, particles[i].n0), particles[i].position);
		tempUArray.push_back(tempU);

		vec3 tempPos = particles[i].position + tempU * mpsTool->GetDeltaT();;
		tempPosArray.push_back(tempPos);

		//用现时刻的位置更新一下粒子的密度
		particles[i].tho = mpsTool->DensityN(posArray, i);
	}
	//用临时速度和临时位置计算每个粒子对应的右端项
	for (int i = 0; i < particles.size(); i++)
	{
		float resDivergence = mpsTool->ExplicitDivergence(tempUArray, tempPosArray, i, particles[i].n0);
		float resRight = mpsTool->ImplicitLaplacianRight(particles[i].n0, resDivergence, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
		Right.push_back(resRight);

		//因为表面检测与右端项的计算不影响，放在同一个循环中提高效率
		particles[i].SurfaceAdjudge(a, g, l0);
		surfaceJudgeArray.push_back(particles[i].isSurface);
		n0Array.push_back(particles[i].n0);
	}

	//1.2.2 解一个泊松方程
	vector<double> resP = mpsTool->ImplicitCalculateP(posArray, n0Array, surfaceJudgeArray, Right);


	//计算每一个粒子新的速度U
	for (int i = 0; i < particles.size(); i++)
	{
		mat3 C = mpsTool->GetMaterixC(posArray, i, particles[i].n0);
		vec3 v1 = mpsTool->ExplicitGradient(C, resP, posArray, particles[i].n0, i);			//计算每个粒子的显式梯度
		vec3 resLU = mpsTool->ExplicitLaplacian(uArray, posArray, i, particles[i].n0);
		particles[i].speed = mpsTool->CalculateU(resLU, v1, particles[i].speed, particles[i].tho);
	}

	//更新位置和压力
	for (int i = 0; i < particles.size(); i++)
	{
		particles[i].pressure = resP[i];
		particles[i].position += particles[i].speed * mpsTool->GetDeltaT();
	}

	//将计算好的位置点进行建模
	//Modeling();
}

void MPSWaterParticleGroup::Draw()
{

}

void MPSWaterParticleGroup::InitBufferData()
{
}
