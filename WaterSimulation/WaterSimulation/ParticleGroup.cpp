#include "ParticleGroup.h"

void MPSWaterParticleGroup::InitParticles()
{
	//初始化MPS工具
	InitMPSTool();
	ParticlesArrange();
	SetInitialN0();

	//更新邻接点关系
	//UpdateAdjoin(range);
	//设置每一个粒子的初始密度
	//SetInitialN0();


	////此处计算所有粒子的初始压力
	//MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	//vector<double> Right;					//存储每一个粒子的右端项

	//vector<vec3> posArray;					//记录k时刻所有粒子的位置
	//vector<vec3> posArray1;					//记录k时刻除去dummy的粒子的位置
	//vector<vec3> uArray;					//记录k时刻所有粒子的速度
	//vector<vec3> uArray1;					//记录k时刻除去dummy的粒子的速度
	//vector<vec3> tempUArray;				//记录忽略压力项后k+1时刻所有粒子的u*
	//vector<vec3> tempUArray1;				//记录忽略压力项后k+1时刻除去dummy的粒子的u*
	//vector<vec3> tempPosArray;				//记录忽略压力项后k+1时刻所有粒子的r*
	//vector<vec3> tempPosArray1;				//记录忽略压力项后k+1时刻除去dummy的粒子的r*

	////1.2 计算隐式的P
	//vector<bool> surfaceJudgeArray;			//记录粒子的表面判断
	//vector<float> n0Array;					//记录粒子的初始密度

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
	////遍历粒子 获取临时的速度和位置（只有普通粒子需要计算，且计算的时候不需要dummy的相关数据,而wall粒子的速度和位置都不变）
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	if (i < particleNumber)
	//	{
	//		//计算临时速度u* 并更新临时位置r*
	//		vec3 tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(uArray1, posArray1, i, particles[i].n0), particles[i].speed);
	//		tempUArray1.push_back(tempU);

	//		vec3 tempPos = particles[i].position + tempU * mpsTool->GetDeltaT();
	//		tempPosArray1.push_back(tempPos);
	//	}
	//	else
	//	{
	//		tempUArray1.push_back(vec3(0));						//wall粒子速度一直为0
	//		tempPosArray1.push_back(particles[i].position);		//wall粒子位置保持不变
	//	}
	//}
	////为全部粒子临时值集合赋值
	//tempUArray = tempUArray1;
	//tempPosArray = tempPosArray1;
	//for (int i = (particleNumber + particleWallNum); i < particleTotalNum; i++)
	//{
	//	tempUArray.push_back(vec3(0));
	//	tempPosArray.push_back(particles[i].position);
	//}


	////用临时速度和临时位置计算每个粒子对应的右端项
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	//普通粒子和wall粒子的右端项计算需要不同的粒子集合
	//	//初始化的时候用旧的右端项方法
	//	if (i < particleNumber)
	//	{
	//		//普通粒子不需要dummy
	//		float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray1, i));
	//		Right.push_back(resRight);
	//	}
	//	else
	//	{
	//		//wall 粒子用全部粒子集合
	//		float resRight = mpsTool->OldImplicitLaplacianRight(particles[i].n0, particles[i].n0, mpsTool->DensityN(tempPosArray, i));
	//		Right.push_back(resRight);
	//	}
	//	
	//	//因为表面检测与右端项的计算不影响，放在同一个循环中提高效率，初始化的时候不需要再进行一次表面检测
	//	surfaceJudgeArray.push_back(particles[i].isSurface);
	//	n0Array.push_back(particles[i].n0);
	//}

	////1.2.2 解一个泊松方程(此时的计算排除dummy)
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
	vec3 offset = vec3(l0);									//粒子空间分布布局
	int width, height, depth;								//粒子的长宽高排布
	width = height = depth = 4;

	height = 5;

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
			R.particleType = SURFACE;
		else
			R.particleType = FLUID;
		R.index = i;
		particles.push_back(R);
	}
	particleNumber = n;

	//还需要初始化粒子墙以及dummy wall	(可以看成再创建1+2层空心无盖的cube)   ps:要确保容器包围住所有的粒子
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
		R.index = i + particleNumber;				//加上之前已经有的索引偏移量
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
		R.index = i + particleNumber;				//加上之前已经有的索引偏移量
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
	//第二层
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
		R.index = i + particleNumber + particleWallNum;				//加上之前已经有的索引偏移量
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

// 计算对应的3个n0值，计算方法为假设出一个被完全包围在中间的粒子（粒子影响范围内充满粒子），用这个粒子的密度作为初始密度
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



	//vector<vec3> posArray;								//所有粒子的位置集合，主要是wall粒子需要
	//vector<vec3> posArray1;								//去除dummy粒子的集合，用于普通粒子
	//for (int i = 0; i < particles.size(); i++)
	//{
	//	posArray.push_back(particles[i].position);
	//	if (i < (particleNumber + particleWallNum))
	//		posArray1.push_back(particles[i].position);
	//}

	////只需要计算普通粒子和墙粒子的密度
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	if (i < particleNumber)			//普通粒子
	//		particles[i].n0 = tool->DensityN(posArray1, i);
	//	else							//墙粒子
	//		particles[i].n0 = tool->DensityN(posArray, i);
	//}

}

//void MPSWaterParticleGroup::UpdateAdjoin(float range)
//{
//	int num1 = particleNumber + particleWallNum;
//	int num2 = particles.size();
//	//邻接粒子中普通粒子不需要考虑dummy，wall粒子都要考虑
//	for (int i = 0; i < num1; i++)
//	{
//		particles[i].adjoinParticleIndex.clear();
//		if (i < particleNumber)
//		{
//			//普通粒子的处理
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
//			//wall粒子的处理
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
	//计算显示步骤（忽略压力项的第一次速度坐标校正）

	vector<vec3> initUArrayWallAndFluid;			// 时间步n的速度集合
	vector<vec3> initPArrayWallAndFluid;			// 时间步n的位置集合
	vector<vec3> middleUArray;						//
	vector<vec3> middlePArray;						//

	vector<vec3> middlePAllArray;					// 所有粒子的中间值坐标（用来计算Wall的中间密度）

	vector<double> Right;					//存储每一个粒子的右端项
	vector<bool> ignoreArray;				// 压力泊松方程中需要忽略的项（surface）

	for (int i = 0; i < particleNumber + particleWallNum; i++)
	{
		initUArrayWallAndFluid.push_back(particles[i].speed);
		initPArrayWallAndFluid.push_back(particles[i].position);
	}

	// 计算每个流体粒子的中间速度
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

	// 将dummy补齐
	for (int i = particleNumber + particleWallNum; i < particleTotalNum; i++)
		middlePAllArray.push_back(particles[i].position);

	// Fluid的中间粒子数密度
	for (int i = 0; i < particleNumber; i++)
	{
		float tempN = mpsTool->DensityN(middlePArray, i);
		particles[i].middleN = tempN;
		//cout << n0ForNumberDensity << " " << tempN << endl;
		// 表面检测
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
	// Wall的中间粒子数密度
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


	// Wall和Fluid粒子的右端项（其中要排除对surface的计算）
	for (int i = 0; i < particleNumber + particleWallNum; i++)
	{
		float resRight = 0;
		//if (!ignoreArray[i])
		//{

		//}
		resRight = mpsTool->OldImplicitLaplacianRight(FLUID_DENSITY, n0ForNumberDensity, particles[i].middleN);
		Right.push_back(resRight);
	}



	// 解一个泊松方程
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
	//每一帧的主要算法流程遍历每一个粒子

	//1.计算每个粒子的p
	//1.1计算每个粒子的右端项
	// u*的散度
	// u的拉普拉斯结果
	// u*的计算
	// n*的计算
	//MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	//vector<double> Right;					//存储每一个粒子的右端项

	//vector<vec3> posArray;					//记录k时刻所有粒子的位置
	//vector<vec3> posArray1;					//记录k时刻除去dummy的粒子的位置
	//vector<vec3> uArray;					//记录k时刻所有粒子的速度
	//vector<vec3> uArray1;					//记录k时刻除去dummy的粒子的速度
	//vector<vec3> tempUArray;				//记录忽略压力项后k+1时刻所有粒子的u*
	//vector<vec3> tempUArray1;				//记录忽略压力项后k+1时刻除去dummy的粒子的u*
	//vector<vec3> tempPosArray;				//记录忽略压力项后k+1时刻所有粒子的r*
	//vector<vec3> tempPosArray1;				//记录忽略压力项后k+1时刻除去dummy的粒子的r*

	////1.2 计算隐式的P
	//vector<bool> surfaceJudgeArray;			//记录粒子的表面判断
	//vector<float> n0Array;					//记录粒子的初始密度

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
	////遍历粒子 获取临时的速度和位置（只有普通粒子需要计算，且计算的时候不需要dummy的相关数据,而wall粒子的速度和位置都不变）
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	//计算临时值的时候顺便计算一下现时刻的粒子密度
	//	if (i < particleNumber)
	//	{
	//		//计算临时速度u* 并更新临时位置r*
	//		vec3 tempU = mpsTool->TempU(mpsTool->ExplicitLaplacian(uArray1, posArray1, i, particles[i].n0), particles[i].speed);
	//		tempUArray1.push_back(tempU);

	//		vec3 tempPos = particles[i].position + tempU * mpsTool->GetDeltaT();
	//		tempPosArray1.push_back(tempPos);

	//		particles[i].tho = mpsTool->DensityN(posArray1, i);
	//	}
	//	else
	//	{
	//		tempUArray1.push_back(vec3(0));						//wall粒子速度一直为0
	//		tempPosArray1.push_back(particles[i].position);		//wall粒子位置保持不变
	//		particles[i].tho = mpsTool->DensityN(posArray, i);
	//	}
	//}
	////为全部粒子临时值集合赋值
	//tempUArray = tempUArray1;
	//tempPosArray = tempPosArray1;
	//for (int i = (particleNumber + particleWallNum); i < particleTotalNum; i++)
	//{
	//	tempUArray.push_back(vec3(0));
	//	tempPosArray.push_back(particles[i].position);
	//}

	////vector<float> tempNArray;
	////用临时速度和临时位置计算每个粒子对应的右端项
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	float tempN = 0;
	//	//普通粒子和wall粒子的右端项计算需要不同的粒子集合
	//	if (i < particleNumber)
	//	{
	//		//普通粒子不需要dummy
	//		float resDivergence = mpsTool->ExplicitDivergence(tempUArray1, tempPosArray1, i, particles[i].n0);
	//		tempN = mpsTool->DensityN(tempPosArray1, i);
	//		//tempNArray.push_back(tempN);
	//		float resRight = mpsTool->ImplicitLaplacianRight(particles[i].n0, resDivergence, particles[i].n0, tempN);
	//		Right.push_back(resRight);
	//	}
	//	else
	//	{
	//		//wall 粒子用全部粒子集合
	//		float resDivergence = mpsTool->ExplicitDivergence(tempUArray, tempPosArray, i, particles[i].n0);
	//		tempN = mpsTool->DensityN(tempPosArray, i);
	//		//tempNArray.push_back(tempN);
	//		float resRight = mpsTool->ImplicitLaplacianRight(particles[i].n0, resDivergence, particles[i].n0, tempN);
	//		Right.push_back(resRight);
	//	}

	//	//因为表面检测与右端项的计算不影响，放在同一个循环中提高效率
	//	//检测的时候wall也检测，因为有dummy所以不会被误认为是表面
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

	////1.2.2 解一个泊松方程(此时的计算排除dummy)
	//vector<double> resP = mpsTool->ImplicitCalculateP(posArray1, n0Array, surfaceJudgeArray, Right);

	////计算每一个粒子新的速度U（只需要计算普通粒子就可以，但是wall粒子需要参与各公式）
	//for (int i = 0; i < particleNumber; i++)
	//{
	//	mat3 C = mpsTool->GetMaterixC(posArray1, i, particles[i].n0);
	//	vec3 v1 = mpsTool->ExplicitGradient(C, resP, posArray1, particles[i].n0, i);			//计算每个粒子的显式梯度
	//	vec3 resLU = mpsTool->ExplicitLaplacian(uArray1, posArray1, i, particles[i].n0);
	//	particles[i].speed = mpsTool->CalculateU(resLU, v1, particles[i].speed, particles[i].tho);
	//}

	////更新位置和压力（只有普通粒子的位置和压力有变化，wall粒子只有压力会变）
	//for (int i = 0; i < (particleNumber + particleWallNum); i++)
	//{
	//	if (particles[i].particleType == SURFACE)
	//		particles[i].pressure = 0;
	//	else
	//		particles[i].pressure = resP[i];
	//	particles[i].position += particles[i].speed * mpsTool->GetDeltaT();
	//}

	//将计算好的位置点进行建模
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
