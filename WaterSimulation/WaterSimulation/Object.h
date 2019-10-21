#pragma once
#include<string>
#include<memory>
#include<OpenMesh/Core/IO/MeshIO.hh>
#include<OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
using namespace std;
#include"Renderer.h"
typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

//基类Object
class Object
{
protected:
	string name;									//object名称
	Transform transformation;						//和空间位置有关的transform组件
	ShaderData* shaderData;							//每一个物体的渲染数据，此处为抽象基类，使用不同渲染器时初始化为相应子类
	Renderer* renderer;								//只是一个指针，不同的渲染器都是单例,不同的物体初始化时只需要将此指针赋值就行
protected:
	//void UpdateMatrix() { shaderData->UpdateMatrix(transformation); }
public:
	//Get
	string GetName() { return name; }
	Transform& GetTransform() { return transformation; }
	ShaderData* GetShaderData() { return shaderData; }
	//Set
	void SetName(string _name) { name = _name; }
	void SetRenderer(RENDERERTYPE type);			//设置渲染器并生成对应的shaderData

	virtual void InitBufferData() = 0;
	virtual void Update(float dt) = 0;
	virtual void Draw() = 0;
};

class IGetVertexDataArray
{
public:
	virtual void GetVertexDataArray(vector<float>& data) = 0;
};

class MeshObject:public Object,public IGetVertexDataArray
{
private:
	Mesh mesh;
private:
	void GetVertexDataArray(vector<float>& data) override;
public:
	MeshObject()
	{
		shaderData = NULL;
		renderer = NULL;
	}
	~MeshObject()
	{
		delete shaderData;
	}



	void readObjFile(string fileName);
	void InitBox(float width, float height, float depth);
	void InitSphere(float radius, int slice, int stack);
	void InitGrid(float radius, int slice, int stack);

	void InitBufferData()override;
	void Update(float dt)override;
	void Draw()override;
};



