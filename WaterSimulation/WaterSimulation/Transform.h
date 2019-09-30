//object���������������
#pragma once
#include<glm.hpp>
#include<vector>
using namespace glm;

class Transform
{
public:
	//�������������
	vec3 position;					//�����λ��
	vec3 scaler;					//����ķŴ�ϵ��
	vec3 rotation;					//�������ת
	//float RotateAngle;				//������ת�ĽǶ�

private:


public:
	//���캯��
	Transform(vec3 _pos = vec3(0), vec3 _scaler = vec3(1.0), vec3 _rotation = vec3(0)) :position(_pos), scaler(_scaler), rotation(_rotation)
	{
	}

	void MoveByDir(vec3 dir, float distant);
	void MoveByVector(vec3 displacement);
	void SetPosition(vec3 _position);

	void RotateByAxis(vec3 axis, float angle);
	void SetRotation(vec3 _rotation);
	
	void SetScaler(vec3 _scaler);
};