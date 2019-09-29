#include "Transform.h"

void Transform::MoveByDir(vec3 dir, float distant)
{
	position += distant * dir;
}

void Transform::MoveByVector(vec3 displacement)
{
	position += displacement;
}

void Transform::SetPosition(vec3 _position)
{
	position = _position;
}

void Transform::RotateByAxis(vec3 axis, float angle)
{
	rotation += angle * axis;
}

void Transform::SetRotation(vec3 _rotation)
{
	rotation = _rotation;
}

void Transform::SetScaler(vec3 _scaler)
{
	scaler = _scaler;
}
