#pragma once
#include<string>
#include<memory>
using namespace std;
#include"Renderer.h"
//����Object
class Object
{
protected:
	string name;							//object����
	Transform transformation;				//�Ϳռ�λ���йص�transform���
	ShaderData* shaderData;					//ÿһ���������Ⱦ���ݣ��˴�Ϊ������࣬��ͬ�����ʼ��ʱ��ֵ��Ӧ����
	Renderer* renderer;						//ֻ��һ��ָ�룬��ͬ����Ⱦ�����ǵ���,��ͬ�������ʼ��ʱֻ��Ҫ����ָ�븳ֵ����
public:
	//Get
	string GetName() { return name; }
	Transform& GetTransform() { return transformation; }
	//Set
	void SetName(string _name) { name = _name; }


	virtual void Update() = 0;
	virtual void Draw() = 0;
};