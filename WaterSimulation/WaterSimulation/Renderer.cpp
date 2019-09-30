#include "Renderer.h"

void Renderer::InitProgram(string vt, string ft)
{
	shaderProgram.SetShader(vt.c_str(), ft.c_str());
}

void Renderer::SetTexture(GLuint& texId, int num, GLenum texNum, string samplerName, ShaderProgram& p)
{
	GLuint texLocation;
	glActiveTexture(texNum);							//��������Ԫ(����λ��)��
	glBindTexture(GL_TEXTURE_2D, texId);				//���������󶨵���ǰ���������Ԫ��
	//������ָ����������Ӧ�ĸ�����Ԫ
	texLocation = glGetUniformLocation(p.p, samplerName.c_str());	//��ȡ��������location
	glUniform1i(texLocation, num);									//ָ����������Ӧ��ǰ�󶨵�����Ԫ0
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
	//ShaderData* test = new UE4ShaderData();
	UE4ShaderData* data = dynamic_cast<UE4ShaderData*>(shaderData);
	if (data != NULL)
	{
		glBindVertexArray(data->VAO);				//��ǰ�����úõ�VAO
								//��������任����
		SetUniform("worldViewProj", data->worldViewProj, shaderProgram);
		SetUniform("world", data->world, shaderProgram);
		SetUniform("worldInvTranspose", data->worldInvTranspose, shaderProgram);

		//SetUniform("depthBiasMVP", depthBiasMVP, p);


		//���ݲ����϶������ѡ�񣬽���Ҫ����������ɫ��
		//�Ƚ��Ƿ�ʹ��������shader
		SetUniform("useTexture", data->bUseTexture, shaderProgram);
		if (data->bUseTexture)
		{
			//����������ͼ
			SetUniform("useAlbedo", data->bAlbedo, shaderProgram);
			if (data->bAlbedo)
			{
				SetTexture(data->tAlbedo, 0, GL_TEXTURE0, "albedoMap", shaderProgram);

			}

			//������ͼ
			SetUniform("useNormal", data->bNormal, shaderProgram);
			if (data->bNormal)
			{
				SetTexture(data->tNormal, 1, GL_TEXTURE1, "normalMap", shaderProgram);

			}

			//��������ͼ
			SetUniform("useMetallic", data->bMetallic, shaderProgram);
			if (data->bMetallic)
			{
				SetTexture(data->tMetallic, 2, GL_TEXTURE2, "metallicMap", shaderProgram);

			}
			else
			{
				//�˴���ʱֱ�ӽ�û�н�����ͼ�Ľ����������ִ���shader
				SetUniform("metallicN", 0.2f, shaderProgram);
			}

			//�ֲ���ͼ
			SetUniform("useRoughness", data->bRoughness, shaderProgram);
			if (data->bRoughness)
			{
				SetTexture(data->tRoughness, 3, GL_TEXTURE3, "roughnessMap", shaderProgram);

			}

			//������ao��ͼ
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
