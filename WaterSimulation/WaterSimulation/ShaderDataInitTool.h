#pragma once
#include<GL/glew.h>
#include<vector>
#include<string>
#include<SOIL.h>
using namespace std;

class ShaderDataInitTool
{
private:
	static ShaderDataInitTool* instance;
	ShaderDataInitTool()
	{
	}
public:
	static ShaderDataInitTool* GetShaderDataInitTool()
	{
		if (instance == NULL)
		{
			instance = new ShaderDataInitTool();
		}
		return instance;
	}
	void SetVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexPos);
	void SetVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexPos, vector<float>& vertexNormal);
	void SetVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexPos, vector<float>& vertexNormal, vector<float>& vertexTex);

	void SetTexture(GLuint& texID, string texPath);
};