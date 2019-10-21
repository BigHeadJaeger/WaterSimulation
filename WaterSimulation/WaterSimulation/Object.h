#pragma once
#include<string>
#include<memory>
#include<OpenMesh/Core/IO/MeshIO.hh>
#include<OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
using namespace std;
#include"Renderer.h"
typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

//����Object
class Object
{
protected:
	string name;									//object����
	Transform transformation;						//�Ϳռ�λ���йص�transform���
	ShaderData* shaderData;							//ÿһ���������Ⱦ���ݣ��˴�Ϊ������࣬ʹ�ò�ͬ��Ⱦ��ʱ��ʼ��Ϊ��Ӧ����
	Renderer* renderer;								//ֻ��һ��ָ�룬��ͬ����Ⱦ�����ǵ���,��ͬ�������ʼ��ʱֻ��Ҫ����ָ�븳ֵ����
protected:
	//void UpdateMatrix() { shaderData->UpdateMatrix(transformation); }
public:
	//Get
	string GetName() { return name; }
	Transform& GetTransform() { return transformation; }
	ShaderData* GetShaderData() { return shaderData; }
	//Set
	void SetName(string _name) { name = _name; }
	void SetRenderer(RENDERERTYPE type);			//������Ⱦ�������ɶ�Ӧ��shaderData

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



