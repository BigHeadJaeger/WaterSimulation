#pragma once
#include<glm.hpp>
#include<vector>
using namespace std;
using namespace glm;

class Particle
{
public:
	vec3 position;			//λ��
	vec3 speed;				//���ӵ��ٶ�����
	float pressure;			//�ܵ���ѹ��
	float life;				//����ʱ��
	float size;				//����ʱָ�뾶

	int index;				//��¼���������ֵ
	vector<int> adjoinParticleIndex;			//��¼�ڽ����ӵ�����

public:
	Particle()
	{
		position = vec3(0);
		speed = vec3(0);
		pressure = 0;
		life = 0;
		size = 0;
		index = -1;
	}
	Particle(vec3 _pos, vec3 _speed, float _pressure, float _size = 1, float _life = 0)
	{
		position = _pos;
		speed = _speed;
		pressure = _pressure;
		life = _life;
		size = _size;
		index = -1;
	}

	void UpdateAdjoin(vector<Particle>& particles, float range);
};