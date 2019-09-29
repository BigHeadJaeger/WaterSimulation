#pragma once
#include<string>
#include<memory>
using namespace std;
#include"ShaderData.h"

//基类Object
class Object
{
protected:
	string name;							//object名称
	Transform transformation;				//和空间位置有关的transform组件
public:
	//Get
	string GetName() { return name; }
	Transform& GetTransform() { return transformation; }
	//Set
	void SetName(string _name) { name = _name; }

	virtual void Update() = 0;
	virtual void Draw() = 0;
};