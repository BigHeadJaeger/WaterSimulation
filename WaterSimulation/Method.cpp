#include "stdafx.h"
#include "Method.h"
#include<omp.h>

CMethod::CMethod()			//һЩ�����С�Ŀռ����
{
	m_Plane = new PlaneNode[G_pScene->m_nVertices];
	m_nOutPoints = G_pScene->m_Boundary.size();
	m_nInteralPoints = G_pScene->m_nVertices - m_nOutPoints;

	//��������Ϊ�ڲ������������Ϊȫ��������Ķ�ά����
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

void CMethod::Initialize()				//���߽��ָ����ƽ���4������,�洢��λ����m_plane��n+1��N��λ�ô�
{
	int sideLength = 1;								//�趨ƽ��ı߳�
	int n = m_nOutPoints / 4;						//���߽��ֳ��ķ�
	double m = sideLength / (double)n;				//��һ������һ������ռ�ľ���

	int i;
	int temp;
	//��ָ���ױ�
	for (i = 0,temp=0; temp < n; i++,temp++)
	{
		m_Plane[m_nInteralPoints + i].index = G_pScene->m_Boundary[i];
		m_Plane[m_nInteralPoints + i].position.x = temp*m;
		m_Plane[m_nInteralPoints + i].position.y = 0;
		m_Plane[m_nInteralPoints + i].position.z = 0;
	}

	//�����½��ϵĵ㿪ʼָ���ұ�
	for (i, temp = 0; temp < n; i++, temp++)
	{
		m_Plane[m_nInteralPoints + i].index = G_pScene->m_Boundary[i];
		m_Plane[m_nInteralPoints + i].position.x = sideLength;
		m_Plane[m_nInteralPoints + i].position.y = temp*m;
		m_Plane[m_nInteralPoints + i].position.z = 0;
	}

	//�����Ͻ��ϵĵ㿪ʼָ���ϱ�
	for (i, temp = n; temp > 0; i++, temp--)
	{
		m_Plane[m_nInteralPoints + i].index = G_pScene->m_Boundary[i];
		m_Plane[m_nInteralPoints + i].position.x = temp*m;
		m_Plane[m_nInteralPoints + i].position.y = sideLength;
		m_Plane[m_nInteralPoints + i].position.z = 0;

	}

	int lastN = m_nOutPoints - 3 * n;		//��Ϊǰ��ֳ��ķݵ�ʱ��һ���������ģ��������һ���߰�ʣ��ĵ�ȫ��ָ����
	double lastM = 1 / (double)lastN;

	//�����Ͻ��ϵĵ㿪ʼָ�����
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

	vector<double> b;			//�����ұߵ�ֵ
	vector<double> x;			//��ļ���
	int nRhs = flag;			//b�����x������ÿһ��������������Ԫ�صĸ�������ʱ��2

								//�ȱ��� �ڲ�������� �� �ڲ�������� ��Namdar��������mkl�����ϵ��
	for (int i = 0; i < m_nInteralPoints; i++)
	{
		bool isFirst = true;			//����ǵ�ǰ�еĵ�һ������Ԫ��
		for (int j = 0; j < m_nInteralPoints; j++)
		{
			if ((i != j) && (Namdar[i][j] != 0))
			{
				a.push_back(-Namdar[i][j]);				//��ӷ���ֵ��a��
				ja.push_back(j+1);						//��ӵ�ǰ����ֵ���кţ���Ϊ�к��Ǵ�1��ʼ������j+1��
				if (isFirst)
				{
					ia.push_back(a.size());				//���뵱ǰ�е�һ��Ԫ����a�е���������Ϊ������1��ʼ�����Ծ͵��ڵ�ǰa��size
					isFirst = false;
				}
			}
			else if (i == j)							//
			{
				a.push_back(1);							//���1��a��
				ja.push_back(j + 1);					//��ӵ�ǰ����ֵ���к�
				if (isFirst)
				{
					ia.push_back(a.size());				//���뵱ǰ�е�һ��Ԫ����a�е���������Ϊ������1��ʼ�����Ծ͵��ڵ�ǰa��size
					isFirst = false;
				}
			}
		}
	}
	ia.push_back(a.size() + 1);							//��ia������һ�� ������ֵ����+1��


	b.resize((ia.size() - 1)*nRhs);
														//���췽���ұ���
	for (int i = 0; i < m_nInteralPoints; i++)
	{
		double sum_u = 0;
		double sum_v = 0;
		for (int j = m_nInteralPoints; j < G_pScene->m_nVertices; j++)		//��ʱҪ�ӵ�n+1��N��ʼ
		{

			sum_u += Namdar[i][j] * m_Plane[j].position.x;
			sum_v += Namdar[i][j] * m_Plane[j].position.y;
		}
		b[i] = sum_u;
		b[i + m_nInteralPoints] = sum_v;			//�Ҷ���Ĵ洢��ʽ���ȴ洢ÿһ��������ĵ�һ��ֵ���ٴ洢�ڶ���ֵ
	}

	x.resize((ia.size() - 1)*nRhs);					//resize x����  ia.size() - 1���Դ�������
	SolveEquation(a, ia, ja, b, x, nRhs);			//�ⷽ��

	//�õ���x�⼯����һά�洢�ģ�����Ҫת������Ӧ��m_Plane�е�ǰm_nInteralPoints���ڲ����position�ϣ�����ָ������
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
	int n = ia.size() - 1;				//�õ�δ֪��Ԫ�ظ���
	const int nnz = ia.back() - 1;		//�õ�����Ԫ�ظ���
	if (ja.size() != nnz || a.size() != nnz || b.size()<nRhs*n || x.size()<nRhs*n)
		return false;

	int mtype = 11;							// ��Ϊ�ⲻ�ԳƵķ�����

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

bool Compare(const PlaneNode& p1, const PlaneNode& p2)		//�����Ƚ������������˭�Ƚ�С
{
	return p1.index < p2.index;
}

void CMethod::Replace()
{
	sort(m_Plane, m_Plane + G_pScene->m_nVertices, Compare);		//�Ƚ�����õ�ƽ�水������С������Ϊplaneǰ�벿�����ڵ㣬��������㣬��m_pVertices�е�˳��һ�£�
	for (int i = 0; i < G_pScene->m_nVertices; i++)
	{
		G_pScene->m_pVertices[i].m_vPosition = m_Plane[i].position;		//������õ�ƽ���ϵĵ�������滻��ԭ���Ķ��㼯����
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


