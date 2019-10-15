#include "stdafx.h"
#include "Method.h"
#include<omp.h>

CMethod::CMethod()			//一些数组大小的空间分配
{
	m_Plane = new PlaneNode[G_pScene->m_nVertices];
	m_nOutPoints = G_pScene->m_Boundary.size();
	m_nInteralPoints = G_pScene->m_nVertices - m_nOutPoints;

	//创建行数为内部点个数，列数为全部点个数的二维矩阵
	Namdar = new double*[G_pScene->m_nVertices];
	for (int i = 0; i < G_pScene->m_nVertices; i++)
	{
		Namdar[i] = new double[G_pScene->m_nVertices];
	}	
}

CMethod::~CMethod()
{
	if (m_Plane!=nullptr)
		delete[] m_Plane;

	if (Namdar != nullptr)
	{
		for (int i = 0; i < m_nInteralPoints; i++)
		{
			delete[] Namdar[i];
		}
		delete[] Namdar;
	}
}

void CMethod::Initialize()				//将边界点指定到平面的4条边上,存储的位置是m_plane的n+1到N的位置处
{
	int sideLength = 1;								//设定平面的边长
	int n = m_nOutPoints / 4;						//将边界点分成四份
	double m = sideLength / (double)n;				//在一条边上一个点所占的距离

	int i;
	int temp;
	//先指定底边
	for (i = 0,temp=0; temp < n; i++,temp++)
	{
		m_Plane[m_nInteralPoints + i].index = G_pScene->m_Boundary[i];
		m_Plane[m_nInteralPoints + i].position.x = temp*m;
		m_Plane[m_nInteralPoints + i].position.y = 0;
		m_Plane[m_nInteralPoints + i].position.z = 0;
	}

	//从右下角上的点开始指定右边
	for (i, temp = 0; temp < n; i++, temp++)
	{
		m_Plane[m_nInteralPoints + i].index = G_pScene->m_Boundary[i];
		m_Plane[m_nInteralPoints + i].position.x = sideLength;
		m_Plane[m_nInteralPoints + i].position.y = temp*m;
		m_Plane[m_nInteralPoints + i].position.z = 0;
	}

	//从右上角上的点开始指定上边
	for (i, temp = n; temp > 0; i++, temp--)
	{
		m_Plane[m_nInteralPoints + i].index = G_pScene->m_Boundary[i];
		m_Plane[m_nInteralPoints + i].position.x = temp*m;
		m_Plane[m_nInteralPoints + i].position.y = sideLength;
		m_Plane[m_nInteralPoints + i].position.z = 0;

	}

	int lastN = m_nOutPoints - 3 * n;		//因为前面分成四份的时候不一定是整除的，所以最后一条边把剩余的点全部指定上
	double lastM = 1 / (double)lastN;

	//从左上角上的点开始指定左边
	for (i, temp = lastN; temp > 0; i++, temp--)
	{
		m_Plane[m_nInteralPoints + i].index = G_pScene->m_Boundary[i];
		m_Plane[m_nInteralPoints + i].position.x = 0;
		m_Plane[m_nInteralPoints + i].position.y = temp*lastM;
		m_Plane[m_nInteralPoints + i].position.z = 0;
	}
}

void CMethod::ConstructEquation(int flag)
{
	vector<double> a;			//Values
	vector<int> ia;				//rowIndex
	vector<int> ja;				//columns

	vector<double> b;			//方程右边的值
	vector<double> x;			//解的集合
	int nRhs = flag;			//b数组和x数组中每一个向量所包含的元素的个数，此时是2

								//先遍历 内部点个数行 和 内部点个数列 的Namdar，来构建mkl所需的系数
	for (int i = 0; i < m_nInteralPoints; i++)
	{
		bool isFirst = true;			//标记是当前行的第一个非零元素
		for (int j = 0; j < m_nInteralPoints; j++)
		{
			if ((i != j) && (Namdar[i][j] != 0))
			{
				a.push_back(-Namdar[i][j]);				//添加非零值在a中
				ja.push_back(j+1);						//添加当前非零值的列号（因为列号是从1开始的所以j+1）
				if (isFirst)
				{
					ia.push_back(a.size());				//放入当前行第一个元素在a中的索引，因为索引从1开始，所以就等于当前a的size
					isFirst = false;
				}
			}
			else if (i == j)							//
			{
				a.push_back(1);							//添加1在a中
				ja.push_back(j + 1);					//添加当前非零值的列号
				if (isFirst)
				{
					ia.push_back(a.size());				//放入当前行第一个元素在a中的索引，因为索引从1开始，所以就等于当前a的size
					isFirst = false;
				}
			}
		}
	}
	ia.push_back(a.size() + 1);							//在ia的最后加一个 （非零值个数+1）


	b.resize((ia.size() - 1)*nRhs);
														//构造方程右边数
	for (int i = 0; i < m_nInteralPoints; i++)
	{
		double sum_u = 0;
		double sum_v = 0;
		for (int j = m_nInteralPoints; j < G_pScene->m_nVertices; j++)		//此时要从第n+1到N开始
		{

			sum_u += Namdar[i][j] * m_Plane[j].position.x;
			sum_v += Namdar[i][j] * m_Plane[j].position.y;
		}
		b[i] = sum_u;
		b[i + m_nInteralPoints] = sum_v;			//右端项的存储格式是先存储每一个点坐标的第一个值，再存储第二个值
	}

	x.resize((ia.size() - 1)*nRhs);					//resize x集合  ia.size() - 1可以代表行数
	SolveEquation(a, ia, ja, b, x, nRhs);			//解方程

	//得到的x解集合是一维存储的，现在要转换到相应的m_Plane中的前m_nInteralPoints个内部点的position上，并且指定索引
	int k = 0;
	for (int i = 0; i < m_nInteralPoints; i++)
	{
		m_Plane[i].index = newIndex[i];
		m_Plane[i].position.x = x[i];
		m_Plane[i].position.y = x[i+ m_nInteralPoints];
		m_Plane[i].position.z = 0;
	}
}

bool CMethod::SolveEquation(vector<double> &a, vector<int> &ia, vector<int> &ja,
	vector<double> &b, vector<double> &x, int nRhs)
{
	int n = ia.size() - 1;				//得到未知量元素个数
	const int nnz = ia.back() - 1;		//得到非零元素个数
	if (ja.size() != nnz || a.size() != nnz || b.size()<nRhs*n || x.size()<nRhs*n)
		return false;

	int mtype = 11;							// 设为解不对称的方程组

	std::vector<void*> pt(64, NULL);		// Internal solver memory pointer
	std::vector<int> iparm(64, 0);			// Pardiso control parameters
	int maxfct, mnum, phase, error, msglvl = 0;
	/* Auxiliary variables. */
	double ddum;			// Double dummy
	int idum;				// Integer dummy
	iparm[0] = 1;			// No solver default
	iparm[1] = 2;			// Fill-in reordering from METIS */
	iparm[2] = omp_get_max_threads();			// omp_get_max_threads();	/* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[7] = 2;			// Max numbers of iterative refinement steps
	iparm[9] = 13;			// Perturb the pivot elements with 1E-13
	iparm[10] = 1;			// Use nonsymmetric permutation and scaling MPS
	iparm[17] = -1;			// Output: Number of nonzeros in the factor LU
	iparm[18] = -1;			// Output: Mflops for LU factorization
	iparm[19] = 0;			// Output: Numbers of CG Iterations
	maxfct = 1;				// Maximum number of numerical factorizations
	mnum = 1;				// Which factorization to use
							//	msglvl = 1;				// Print statistical information in file
	error = 0;				// Initialize error flag

							//////////////////////////////////////////////////////////////////////////
							// .. Reordering and Symbolic Factorization. This step also allocates
							// all memory that is necessary for the factorization. */
							//////////////////////////////////////////////////////////////////////////
	phase = 11;
	PARDISO(&pt.front(), &maxfct, &mnum, &mtype, &phase, &n, &a.front(), &ia.front(), &ja.front(),
		&idum, &nRhs, &iparm.front(), &msglvl, &ddum, &ddum, &error);

	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		return false;
	}
	//////////////////////////////////////////////////////////////////////////
	// .. Numerical factorization
	//////////////////////////////////////////////////////////////////////////
	phase = 22;
	PARDISO(&pt.front(), &maxfct, &mnum, &mtype, &phase, &n, &a.front(), &ia.front(), &ja.front(),
		&idum, &nRhs, &iparm.front(), &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		return false;
	}

	//////////////////////////////////////////////////////////////////////////
	// .. Back substitution and iterative refinement
	//////////////////////////////////////////////////////////////////////////
	phase = 33;
	PARDISO(&pt.front(), &maxfct, &mnum, &mtype, &phase, &n, &a.front(), &ia.front(), &ja.front(),
		&idum, &nRhs, &iparm.front(), &msglvl, &b.front(), &x.front(), &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		return false;
	}

	//////////////////////////////////////////////////////////////////////////
	// .. Termination and release of memory
	//////////////////////////////////////////////////////////////////////////
	phase = -1; /* Release internal memory. */
	PARDISO(&pt.front(), &maxfct, &mnum, &mtype, &phase, &n, &ddum, &ia.front(), &ja.front(),
		&idum, &nRhs, &iparm.front(), &msglvl, &ddum, &ddum, &error);

	return true;
}

bool Compare(const PlaneNode& p1, const PlaneNode& p2)		//用来比较两个点的索引谁比较小
{
	return p1.index < p2.index;
}

void CMethod::Replace()
{
	sort(m_Plane, m_Plane + G_pScene->m_nVertices, Compare);		//先将构造好的平面按索引大小排序（因为plane前半部分是内点，后面是外点，与m_pVertices中的顺序不一致）
	for (int i = 0; i < G_pScene->m_nVertices; i++)
	{
		G_pScene->m_pVertices[i].m_vPosition = m_Plane[i].position;		//将构造好的平面上的点的坐标替换到原来的顶点集合中
	}
}

bool CMethod::judge(int v1, int v2)
{
	int i;
	for (i = 0; i < G_pScene->m_pVertices[v1].m_nNeighbor; i++)
	{
		if (G_pScene->m_pVertices[v1].m_npNeighborVertexIndices[i] == v2)
			break;
	}
	if (i == G_pScene->m_pVertices[v1].m_nNeighbor)
		return false;
	else
		return true;
}

bool CMethod::exist(vector<int> Array, int flag)
{
	for (int i = 0; i < Array.size(); i++)
	{
		if (Array[i] == flag)
			return true;
	}
	return false;
}


