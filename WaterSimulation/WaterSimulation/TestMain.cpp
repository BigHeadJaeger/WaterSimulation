//���㷨���̵�α����
#include"Particle.h"
vector<Particle> particles;

float particleNumber = 10;
float l0 = 1;
float range = 0.5;
float viscosity = 1;
float a = 0.75;


int main()
{
	/*************��ʼ���׶�****************/
	//��Mps���̹��ߵ�һЩ��ֵ���г�ʼ��
	MPSToolFun* mpsTool = MPSToolFun::GetMPSTool();
	mpsTool->SetRe(l0);
	mpsTool->SetViscosity(viscosity);
	//���ӵĳ�ʼ�����������Եĳ�ʼ�������뵽�����б���
	for (int i = 0; i < particleNumber; i++)
	{
		//��ʼ����ǰ���Ӳ�push��ȥ
		Particle R;
		R.index = i;
		particles.push_back(R);
	}
	//����������������Գ�ʼ����֮���ٳ�ʼ������֮����Ӱ�������
	for (int i = 0; i < particles.size(); i++)
	{
		//��һ���ڽ����ӵ��������ϳ�ʼ��
		particles[i].UpdateAdjoin(particles, range);
		//��ȡ�ڽӵ��λ�ü���
		vector<vec3> neighborPos;
		neighborPos.push_back(particles[i].position);
		for (int j = 0; j < particles[i].adjoinParticleIndex.size(); j++)
		{
			int neighborIndex = particles[i].adjoinParticleIndex[j];
			neighborPos.push_back(particles[neighborIndex].position);
		}
		//���õ�ǰ��ĳ�ʼ�ܶ�
		particles[i].SetInitialN0(neighborPos, particles[i].index);
	}


	/*************ÿһ֡������****************/
	float deltaT = 1 / 60;
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
		vector<vec3> neighborPos;
		vector<vec3> neighborU;
		neighborPos.push_back(particles[i].position);
		neighborU.push_back(particles[i].speed);
		vec3 tempU = mpsTool->TempU(deltaT, mpsTool->ExplicitLaplacian(neighborU, neighborPos, i, particles[i].n0), particles[i].position);
		Right.push_back(mpsTool->OldImplicitLaplacianRight(particles[i].n0, deltaT, particles[i].n0, mpsTool->DensityN(neighborPos, i)));
	}
	//1.2 ������ʽ��P
	//1.2.1 ���������ӽ��б�����
	for (int i = 0; i < particles.size(); i++)
	{
		vector<vec3> neighborPos;
		neighborPos.push_back(particles[i].position);
		particles[i].SurfaceAdjudge(a, mpsTool->DensityN(neighborPos, i), g, l0);
	}
	//1.2.2 ��һ�����ɷ���



}