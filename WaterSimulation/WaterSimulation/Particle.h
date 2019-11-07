#pragma once
#include"MPSFun.h"

//ֻ������������Եĳ�������
class Particle
{
public:
	vec3 position;
	bool state;
	float life;									//����ʱ��
	float size;									//����ʱָ�뾶
};

class MPSWaterParticle:public Particle
{
public:
	vec3 speed;									//���ӵ��ٶ�����
	float pressure;								//�ܵ���ѹ��
	float n0;									//��ʼ���ܶ�ֵ
	float tho;									//��ǰ���ܶ�

	bool isSurface;								//�����Ƿ�Ϊ����
	bool isWall;								//�����Ƿ�Ϊǽ
	bool isDummy;								//�����Ƿ�Ϊ����ǽ�棨����������ǽ���Ӻͱ������ӣ�

	int index;									//��¼���������ֵ
	vector<int> adjoinParticleIndex;			//��¼�ڽ����ӵ�����

public:
	MPSWaterParticle()
	{
		position = vec3(0);
		speed = vec3(0);
		pressure = 0;
		n0 = 0;
		tho = 0;
		life = 0;
		size = 0;
		index = -1;
		isSurface = false;
		isWall = false;
		isDummy = false;
	}
	MPSWaterParticle(vec3 _pos, vec3 _speed, float _pressure, float _size = 1, float _life = 0)
	{
		position = _pos;
		speed = _speed;
		pressure = _pressure;
		life = _life;
		size = _size;
		index = -1;
	}

	void SurfaceAdjudge(float a, float g, float l0);
};