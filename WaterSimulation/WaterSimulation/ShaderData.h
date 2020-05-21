#pragma once
//��OpenGLshader������صĸ��ֱ����ͷ���

#include"ShaderDataInitTool.h"
#include"Transform.h"
#include"PublicStruct.h"
#include"Camera.h"
using namespace std;



class ShaderData
{
public:
	//�������������
	mat4x4 world;					//�������
	mat4x4 worldViewProj;			//��������ת������
	mat4x4 worldInvTranspose;		//������������ת��������ռ�

	//�����VAO��VBO���
	GLuint VAO;
	GLuint VBO;
	GLuint IndexBuffer;

	GLint drawType;					//����buffer�Ļ��Ʒ�ʽ
	GLint drawUnitNumber;			//���Ƶ�Ԫ������
public:
	ShaderData()
	{
		world = mat4(0);
		worldViewProj = mat4(0);
		worldInvTranspose = mat4(0);
		drawType = GL_TRIANGLES;
	}

	void SetDrawType(GLenum type) { drawType = type; }

	void UpdateMatrix(Transform& t);
	void InitVertexBuffer(vector<float>& vertexData, bool providedNormal, bool providedTex, GLenum usage = GL_STATIC_DRAW);
	void UpdateVertexBuffer(vector<float>& vertexData, bool providedNormal, bool providedTex);
	virtual void Temp(){}
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
public:
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
	void InitTexture(TEXTURETYPE type, string texPath);
};

class SimpleShaderData :public ShaderData
{
public:
	vec3 color;
public:
	void SetColor(vec3 _color)
	{
		color = _color;
	}
};

class VCShaderData : public ShaderData
{
public:
	void InitVertexBuffer(vector<float>& vertexData, GLenum usage = GL_STATIC_DRAW);
};