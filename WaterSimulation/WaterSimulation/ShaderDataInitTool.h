#pragma once
#include<GL/glew.h>
#include<vector>
#include<string>
#include<SOIL.h>
#include<iostream>
using namespace std;

class ShaderDataInitTool
{
private:
	static ShaderDataInitTool* instance;
	ShaderDataInitTool()
	{
	}

	// ������Ķ�����Ϣ�Ļ�������ɢ��pos normal tex����������
	void VertexDataUnpack(vector<float>& vertexData, vector<float>& vertexPos, vector<float>& vertexNormal, vector<float>& vertexTex);
public:
	static ShaderDataInitTool* GetShaderDataInitTool()
	{
		if (instance == NULL)
		{
			instance = new ShaderDataInitTool();
		}
		return instance;
	}
	void InitVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexData, bool providedNormal, bool providedTex, GLenum usage = GL_STATIC_DRAW);
	void UpdateVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexData, bool providedNormal, bool providedTex);

	void InitTextureWithFile(GLuint& texID, string texPath);


	void InitBufferVC(GLuint& VAO, GLuint& VBO, vector<float>& vertexData, GLenum usage = GL_STATIC_DRAW);
};