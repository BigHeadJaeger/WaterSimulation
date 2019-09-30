#pragma once
#include<glm.hpp>
#include<gtc\matrix_transform.hpp>
#include<gtc\type_ptr.hpp>
using namespace glm;
using namespace std;
#include"Program.h"
#include"ShaderData.h"

class Renderer
{
protected:
	ShaderProgram shaderProgram;
public:
	void InitProgram(string vt, string ft);

	virtual void Render(ShaderData* shaderData) = 0;

	//传texture到shader中
	void SetTexture(GLuint& texId, int num, GLenum texNum, string samplerName, ShaderProgram& p);
	//根据不同类型的值用重载的方式传入shader中
	void SetUniform(string valueName, mat4x4& value, ShaderProgram& p);
	void SetUniform(string valueName, vec4& value, ShaderProgram& p);
	void SetUniform(string valueName, vec3& value, ShaderProgram& p);
	void SetUniform(string valueName, float value, ShaderProgram& p);
};

//不同的渲染器只需要一个，所以都设为单例
class UE4Renderer:public Renderer
{
private:
	static UE4Renderer* instance;
	UE4Renderer(){}
public:
	static UE4Renderer* GetRenderer()
	{
		if (instance == NULL)
		{
			instance = new UE4Renderer();
		}
		return instance;
	}
	void Render(ShaderData* shaderData) override;
};

class MPSRenderer:public Renderer
{
private:
	static MPSRenderer* instance;
	MPSRenderer() {}
public:
	static MPSRenderer* GetRenderer()
	{
		if (instance == NULL)
		{
			instance = new MPSRenderer();
		}
		return instance;
	}
	void Render(ShaderData* shaderData) override;
};