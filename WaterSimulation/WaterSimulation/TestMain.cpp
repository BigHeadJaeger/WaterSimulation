//���㷨���̵�α����
#include"Particle.h"
vector<Particle> particles;

float particleNumber = 10;
float l0 = 1;
float range = 0.5;
float viscosity = 1;

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
	for (int i = 0; i < particles.size(); i++)
	{
		
	}



}