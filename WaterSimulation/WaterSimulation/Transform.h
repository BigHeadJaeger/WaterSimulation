//object的坐标相关属性类
#pragma once
#include<glm.hpp>
#include<vector>
#include"Camera.h"
using namespace glm;

class Transform
{
public:
	//物体的坐标属性
	vec3 position;					//物体的位置
	vec3 scaler;					//物体的放大系数
	vec3 rotation;					//物体的旋转
	//float RotateAngle;				//物体旋转的角度

private:


public:
	//构造函数
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