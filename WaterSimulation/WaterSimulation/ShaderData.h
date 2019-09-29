#pragma once
//��OpenGLshader������صĸ��ֱ����ͷ���

#include"ShaderDataInitTool.h"
#include"Transform.h"
#include"Camera.h"
#include"PublicStruct.h"
using namespace std;



class ShaderData
{
private:


protected:
	//�������������
	mat4x4 world;					//�������
	mat4x4 worldViewProj;			//��������ת������
	mat4x4 worldInvTranspose;		//������������ת��������ռ�

	//�����VAO��VBO���
	GLuint VAO;
	GLuint VertexBuffer;
	GLuint IndexBuffer;


public:
	ShaderData()
	{
		world = mat4(0);
		worldViewProj = mat4(0);
		worldInvTranspose = mat4(0);
	}

	void UpdateMatrix(Transform& t, Camera& camera);
			//����������󲢰�����(��ͼƬת��Ϊ�������ݣ����ݲ�ͬ��ID������Ӧ������)
	void InitVertexBuffer(vector<float>& vertexPos, vector<float>& vertexNormal, vector<float>& vertexTex, bool providedTex);
};

enum TEXTURETYPE
{
	ALBEDO,
	METALLIC,
	ROUGHNESS,
	AO,
	NORMAL,  
};

class UE4ShaderData :public ShaderData
{
protected:
	//�������ͼ��ţ���һ��ȫ����Ҫ��
	bool bUseTexture;
	GLuint tAlbedo;					//������ͼ��������ɫ��
	bool bAlbedo;
	GLuint tMetallic;				//��������ͼ
	bool bMetallic;
	GLuint tRoughness;				//�ֲ���ͼ
	bool bRoughness;
	GLuint tAo;						//������ͼ
	bool bAo;
	GLuint tNormal;					//������ͼ
	bool bNormal;
public:
	void SetTexture(TEXTURETYPE type, string texPath);
};