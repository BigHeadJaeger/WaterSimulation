#include "Renderer.h"

void Renderer::InitProgram(string vt, string ft)
{
	shaderProgram.SetShader(vt.c_str(), ft.c_str());
}

void Renderer::SetTexture(GLuint& texId, int num, GLenum texNum, string samplerName, ShaderProgram& p)
{
	GLuint texLocation;
	glActiveTexture(texNum);							//激活纹理单元(纹理位置)。
	glBindTexture(GL_TEXTURE_2D, texId);				//将纹理对象绑定到当前激活的纹理单元处
	//接下来指定采样器对应哪个纹理单元
	texLocation = glGetUniformLocation(p.p, samplerName.c_str());	//获取采样器的location
	glUniform1i(texLocation, num);									//指定采样器对应当前绑定的纹理单元0
}

void Renderer::SetUniform(string valueName, mat4x4& value, ShaderProgram& p)
{
	GLuint location;
	location = glGetUniformLocation(p.p, valueName.c_str());
	glUniformMatrix4fv(location, 1, GL_FALSE, value_ptr(value));
}

void Renderer::SetUniform(string valueName, vec4& value, ShaderProgram& p)
{
	GLuint location;
	location = glGetUniformLocation(p.p, valueName.c_str());
	glUniform4fv(location, 1, value_ptr(value));
}

void Renderer::SetUniform(string valueName, vec3& value, ShaderProgram& p)
{
	GLuint location;
	location = glGetUniformLocation(p.p, valueName.c_str());
	glUniform3fv(location, 1, value_ptr(value));
}

void Renderer::SetUniform(string valueName, float value, ShaderProgram& p)
{
	GLuint location;
	location = glGetUniformLocation(p.p, valueName.c_str());
	glUniform1f(location, value);
}

void UE4Renderer::Render(ShaderData* shaderData)
{
	glUseProgram(shaderProgram.p);
	//ShaderData* test = new UE4ShaderData();
	UE4ShaderData* data = dynamic_cast<UE4ShaderData*>(shaderData);
	if (data != NULL)
	{
		glBindVertexArray(data->VAO);				//绑定前面设置好的VAO
								//传递坐标变换矩阵
		SetUniform("worldViewProj", data->worldViewProj, shaderProgram);
		SetUniform("world", data->world, shaderProgram);
		SetUniform("worldInvTranspose", data->worldInvTranspose, shaderProgram);

		//SetUniform("depthBiasMVP", depthBiasMVP, p);


		//根据参数上对纹理的选择，将需要的纹理传入着色器
		//先将是否使用纹理传入shader
		SetUniform("useTexture", data->bUseTexture, shaderProgram);
		if (data->bUseTexture)
		{
			//基础反射贴图
			SetUniform("useAlbedo", data->bAlbedo, shaderProgram);
			if (data->bAlbedo)
			{
				SetTexture(data->tAlbedo, 0, GL_TEXTURE0, "albedoMap", shaderProgram);

			}

			//法线贴图
			SetUniform("useNormal", data->bNormal, shaderProgram);
			if (data->bNormal)
			{
				SetTexture(data->tNormal, 1, GL_TEXTURE1, "normalMap", shaderProgram);

			}

			//金属度贴图
			SetUniform("useMetallic", data->bMetallic, shaderProgram);
			if (data->bMetallic)
			{
				SetTexture(data->tMetallic, 2, GL_TEXTURE2, "metallicMap", shaderProgram);

			}
			else
			{
				//此处暂时直接将没有金属贴图的金属度用数字传入shader
				SetUniform("metallicN", 0.2f, shaderProgram);
			}

			//粗糙贴图
			SetUniform("useRoughness", data->bRoughness, shaderProgram);
			if (data->bRoughness)
			{
				SetTexture(data->tRoughness, 3, GL_TEXTURE3, "roughnessMap", shaderProgram);

			}

			//环境光ao贴图
			SetUniform("useAO", data->bAo, shaderProgram);
			if (data->bAo)
			{
				SetTexture(data->tAo, 4, GL_TEXTURE4, "aoMap", shaderProgram);

			}
		}

		glDrawArrays(data->drawType, 0, data->drawUnitNumber);

		glBindVertexArray(0);
	}
}

void SimpleRenderer::Render(ShaderData* shaderData)
{
	glUseProgram(shaderProgram.p);
	SimpleShaderData* data = dynamic_cast<SimpleShaderData*>(shaderData);
	glBindVertexArray(data->VAO);
	//传递坐标变换矩阵
	SetUniform("worldViewProj", data->worldViewProj, shaderProgram);
	SetUniform("world", data->world, shaderProgram);
	SetUniform("worldInvTranspose", data->worldInvTranspose, shaderProgram);

	SetUniform("color", data->color, shaderProgram);
	glDrawArrays(data->drawType, 0, data->drawUnitNumber);
}

void VCRenter::Render(ShaderData* shaderData)
{
	glUseProgram(shaderProgram.p);
	VCShaderData* data = dynamic_cast<VCShaderData*>(shaderData);
	glBindVertexArray(data->VAO);
	//传递坐标变换矩阵
	SetUniform("worldViewProj", data->worldViewProj, shaderProgram);
	SetUniform("world", data->world, shaderProgram);
	SetUniform("worldInvTranspose", data->worldInvTranspose, shaderProgram);
	glDrawArrays(data->drawType, 0, data->drawUnitNumber);
}


UE4Renderer* UE4Renderer::instance = NULL;
SimpleRenderer* SimpleRenderer::instance = NULL;
MPSRenderer* MPSRenderer::instance = NULL;
VCRenter* VCRenter::instance = NULL;

