#include "MPSFun.h"

MPSToolFun* MPSToolFun::mpsTool = NULL;


float MPSToolFun::WeightFun(float dis,float re)
{
	float result = 0;
	if (dis >= 0 && dis < re)
		result = re / dis - 1;
	else if (dis >= re)
		result = 0;
	return result;
}

vec3 MPSToolFun::TempU(vec3 resLU, vec3 uNow)
{
	return deltaT * (viscosity * resLU + vG) + uNow;
}

vec3 MPSToolFun::CalculateU(vec3 resLU, vec3 resGP, vec3 uNow, float tho)
{
	return deltaT * ((-1 / tho) * resGP + viscosity * resLU + vG) + uNow;
}

vector<double> MPSToolFun::ImplicitCalculateP(vector<vec3>& r, vector<float>& n0Array, vector<bool>& isSurface,vector<double> Right)
{
	//ÿһ���е�n0��lambda���ǲ�ͬ��

	vector<double> a;			//Values
	vector<int> ia;				//rowIndex
	vector<int> ja;				//columns

	//vector<double> b;			//�����ұߵ�ֵ
	vector<double> x;			//��ļ���
	int nRhs = 1;			//b�����x������ÿһ��������������Ԫ�صĸ�������ʱ��2

	for (int i = 0; i < n0Array.size(); i++)
	{
		vector<double> coeffArray;			//ϵ������
		coeffArray.resize(n0Array.size(), 0);
		double currentRowCoeff = 0;			//��ǰ�б�δ֪����ϵ����Ҫ����δ֪����ϵ����ͬ����

		bool isFirst = true;				//��ǵ�ǰ���Ƿ��һ������Ԫ��

		double n0 = n0Array[i];				//��ǰ�е�n0
		double lambda = Lambda(r, i);		//��ǰ�е�lambda
		double con = 2 * Ds / (n0 * lambda);	//ÿһ�еĳ���

		for (int j = 0; j < n0Array.size(); j++)
		{	
			if (i != j)
			{
				double w = WeightFun(distance(r[j], r[i]), reForL);
				coeffArray[j] = w * con;
				currentRowCoeff -= coeffArray[j];
			}
		}
		coeffArray[i] = currentRowCoeff;
		//��ǰ�е�ϵ�����������ˣ�����mkl�����ʽ��ӵ���Ӧ������
		for (int k = 0; k < coeffArray.size(); k++)
		{
			if (coeffArray[k] != 0 && !isSurface[k])		//��ǰ���ϵ����Ϊ0�Ҳ��Ǳ����
			{
				a.push_back((coeffArray[k]));
				ja.push_back(k + 1);
				if (isFirst)
				{
					ia.push_back(a.size());
					isFirst = false;
				}
			}
		}
	}

	ia.push_back(a.size() + 1);							//��ia������һ�� (����ֵ���� + 1)

	x.resize((ia.size() - 1));

	bool res = SolveEquation(a, ia, ja, Right, x, 1);
	if(res)
		return x;
	else
	{
		cout << "�������������" << endl;
		return vector<double>();
	}
}

bool MPSToolFun::SolveEquation(vector<double>& a, vector<int>& ia, vector<int>& ja, vector<double>& b, vector<double>& x, int nRhs)
{
	int n = ia.size() - 1;				//�õ�δ֪��Ԫ�ظ���
	const int nnz = ia.back() - 1;		//�õ�����Ԫ�ظ���
	if (ja.size() != nnz || a.size() != nnz || b.size() < nRhs * n || x.size() < nRhs * n)
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

float MPSToolFun::ImplicitLaplacianRight(float rho0, float resDu, float n0, float tempN)
{
	return (1 - gama) * rho0 * (resDu / deltaT) - gama * (rho0 / pow(deltaT, 2)) * ((tempN - n0) / n0);
}

float MPSToolFun::OldImplicitLaplacianRight(float rho0, float n0, float tempN)
{
	return (rho0 / pow(deltaT, 2)) * ((n0 - tempN) / n0);
}

mat3 MPSToolFun::GetMaterixC(vector<vec3>& R, int currentIndex, float n0)
{
	mat3 res(0);
	for (int i = 0; i < R.size(); i++)
	{
		if (i != currentIndex)
		{
			vec3 v1 = (1 / length(R[i] - R[currentIndex])) * (R[i] - R[currentIndex]);
			res += (WeightFun(length(R[i] - R[currentIndex]), reForDG) * outerProduct(v1, v1));
		}
	}
	res *= (1 / n0);
	return res;
}

float MPSToolFun::ExplicitDivergence(vector<vec3>& phi, vector<vec3>& r, int currentIndex, float n0)
{
	float res = 0;
	for (int i = 0; i < phi.size(); i++)
	{
		if (i != currentIndex)
		{
			res += (dot((phi[i] - phi[currentIndex]), (r[i] - r[currentIndex])) / pow(distance(r[i], r[currentIndex]), 2))
				* WeightFun(distance(r[i], r[currentIndex]), reForDG);
		}
	}
	res *= (Ds / n0);
	return res;
}

vec3 MPSToolFun::ExplicitGradient(mat3 C, vector<double>& p, vector<vec3>& r, float n0, int currentIndex)
{
	if (determinant(C) >= 0.05)
	{
		vec3 res(0);
		for (int i = 0; i < p.size(); i++)
		{
			if (i != currentIndex)
			{
				float l = length(r[i] - r[currentIndex]);
				res += (WeightFun(length(l), reForDG) * (((float)p[i] - (float)p[currentIndex]) / l) * ((r[i] - r[currentIndex]) / l));
			}
		}
		res /= n0;
		return inverse(C)* res;
	}
	else
	{
		vec3 res(0);
		for (int i = 0; i < p.size(); i++)
		{
			if (i != currentIndex)
			{
				float l = length(r[i] - r[currentIndex]);
				res += (WeightFun(length(l), reForDG) * (((float)p[i] - (float)p[currentIndex]) / l) * ((r[i] - r[currentIndex]) / l));
			}
		}
		res *= (Ds / n0);
		return res;
	}
}

float MPSToolFun::DensityN(vector<vec3>& r, int currentIndex)
{
	float res = 0;
	for (int i = 0; i < r.size(); i++)
	{
		if (i != currentIndex)
		{
			res += WeightFun(distance(r[i], r[currentIndex]), reForDG);
		}
	}
	return res;
}

vec3 MPSToolFun::NewPosR(vec3 nowPos, vec3 u)
{
	return nowPos + u;
}

float MPSToolFun::Lambda(vector<vec3> r, int currentIndex)
{
	float result = 0;
	float numerator = 0;
	float denominator = 0;
	for (int i = 0; i < r.size(); i++)
	{
		if (i != currentIndex)
		{
			numerator += pow(distance(r[i], r[currentIndex]), 2) * WeightFun(distance(r[i], r[currentIndex]), reForL);
			denominator += WeightFun(distance(r[i], r[currentIndex]), reForL);
		}
	}
	result = numerator / denominator;
	return result;
}
