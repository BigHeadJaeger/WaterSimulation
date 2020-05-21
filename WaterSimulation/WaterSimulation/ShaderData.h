#pragma once
//与OpenGLshader配置相关的各种变量和方法

#include"ShaderDataInitTool.h"
#include"Transform.h"
#include"PublicStruct.h"
#include"Camera.h"
using namespace std;



class ShaderData
{
public:
	//物体的坐标属性
	mat4x4 world;					//世界矩阵
	mat4x4 worldViewProj;			//最终坐标转换矩阵
	mat4x4 worldInvTranspose;		//用来将法向量转换到世界空间

	//物体的VAO、VBO编号
	GLuint VAO;
	GLuint VBO;
	GLuint IndexBuffer;

	GLint drawType;					//顶点buffer的绘制方式
	GLint drawUnitNumber;			//绘制单元的数量
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
	//物体的贴图编号（不一定全都需要）
	bool bUseTexture;
	GLuint tAlbedo;					//反射贴图（基础颜色）
	bool bAlbedo;
	GLuint tMetallic;				//金属度贴图
	bool bMetallic;
	GLuint tRoughness;				//粗糙贴图
	bool bRoughness;
	GLuint tAo;						//环境贴图
	bool bAo;
	GLuint tNormal;					//法线贴图
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