#pragma once
#include<string>
#include<memory>
using namespace std;
#include"Renderer.h"
//基类Object
class Object
{
protected:
	string name;							//object名称
	Transform transformation;				//和空间位置有关的transform组件
	ShaderData* shaderData;					//每一个物体的渲染数据，此处为抽象基类，不同物体初始化时赋值相应子类
	Renderer* renderer;						//只是一个指针，不同的渲染器都是单例,不同的物体初始化时只需要将此指针赋值就行
public:
	//Get
	string GetName() { return name; }
	Transform& GetTransform() { return transformation; }
	//Set
	void SetName(string _name) { name = _name; }


	virtual void Update() = 0;
	virtual void Draw() = 0;
};