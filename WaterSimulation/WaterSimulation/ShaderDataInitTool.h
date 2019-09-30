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
	void InitVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexData, bool providedNormal, bool providedTex);

	void InitTextureWithFile(GLuint& texID, string texPath);
};