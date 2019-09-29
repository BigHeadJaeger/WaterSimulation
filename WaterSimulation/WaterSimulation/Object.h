#pragma once
#include<string>
#include<memory>
using namespace std;
#include"ShaderData.h"

//����Object
class Object
{
protected:
	string name;							//object����
	Transform transformation;				//�Ϳռ�λ���йص�transform���
public:
	//Get
	string GetName() { return name; }
	Transform& GetTransform() { return transformation; }
	//Set
	void SetName(string _name) { name = _name; }

	virtual void Update() = 0;
	virtual void Draw() = 0;
};