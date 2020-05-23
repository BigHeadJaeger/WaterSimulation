#pragma once
#include<map>
using namespace std;
#include"Object.h"
//#include"ParticleGroup.h"
#include"MPSWaterParticleGroup.h"
#include"Interaction.h"

struct DrawMode
{
	bool isLine;
};
class MyScene
{
private:
	//各种场景信息（相机、材质、灯光、各种物体的各种矩阵）
	map<string, Object*> objects;
	//MeshObject cow;
	//灯光
	vec3 lightPos;
	vec3 lightColor;

public:
	//vector<Key> keys;
	map<KEYNAME, Key> keys;

	Mouse mouse;

	DrawMode drawMode;
	//---------------------------------------------------------------------------------

private:

public:
	~MyScene()
	{
		map<string, Object*>::iterator objs_it;
		for (objs_it = objects.begin(); objs_it != objects.end(); objs_it++)
		{
			delete (*objs_it).second;
		}
		objects.clear();
	}

	void Init();			//初始化各种信息
	void InitKeys();

	void Update(float dt);			//需要动画时，计算各种矩阵（暂时不传入shader中）
	void Draw();			//绘制场景
private:
};