#include "Object.h"

void Object::InitBounding(BOUNDINGTYPE type)
{
	switch (type)
	{
	case BOX:
		bounding = make_shared<BoundingBox>();
		bounding->Init(this->meshData);
		break;
	case SPHERE:
		break;
	default:
		break;
	}
}

void Object::Draw(ShaderProgram&p)
{
	ConveyTool* conveyTool = ConveyTool::GetConveyTool();

	glBindVertexArray(shaderData.VAO);				//��ǰ�����úõ�VAO
									//��������任����
	conveyTool->SetUniform("worldViewProj", shaderData.worldViewProj, p);
	conveyTool->SetUniform("world", shaderData.world, p);
	conveyTool->SetUniform("worldInvTranspose", shaderData.worldInvTranspose, p);

	//SetUniform("depthBiasMVP", depthBiasMVP, p);


	//���ݲ����϶������ѡ�񣬽���Ҫ����������ɫ��
	//�Ƚ��Ƿ�ʹ��������shader
	conveyTool->SetUniform("useTexture", shaderData.bUseTexture, p);
	if (shaderData.bUseTexture)
	{
		//������ʹ������������������
		if (meshData.providedTex)
		{
			//����������ͼ
			conveyTool->SetUniform("useAlbedo", shaderData.bAlbedo, p);
			if (shaderData.bAlbedo)
			{
				conveyTool->SetTexture(shaderData.tAlbedo, 0, GL_TEXTURE0, "albedoMap", p);

			}

			//������ͼ
			conveyTool->SetUniform("useNormal", shaderData.bNormal, p);
			if (shaderData.bNormal)
			{
				conveyTool->SetTexture(shaderData.tNormal, 1, GL_TEXTURE1, "normalMap", p);

			}

			//��������ͼ
			conveyTool->SetUniform("useMetallic", shaderData.bMetallic, p);
			if (shaderData.bMetallic)
			{
				conveyTool->SetTexture(shaderData.tMetallic, 2, GL_TEXTURE2, "metallicMap", p);

			}
			else
			{
				//�˴���ʱֱ�ӽ�û�н�����ͼ�Ľ����������ִ���shader
				conveyTool->SetUniform("metallicN", 0.2f, p);
			}

			//�ֲ���ͼ
			conveyTool->SetUniform("useRoughness", shaderData.bRoughness, p);
			if (shaderData.bRoughness)
			{
				conveyTool->SetTexture(shaderData.tRoughness, 3, GL_TEXTURE3, "roughnessMap", p);

			}

			//������ao��ͼ
			conveyTool->SetUniform("useAO", shaderData.bAo, p);
			if (shaderData.bAo)
			{
				conveyTool->SetTexture(shaderData.tAo, 4, GL_TEXTURE4, "aoMap", p);

			}


		}

		//��Ӱ��ͼ����������������ͬһ�����ͣ�����Ҫ����������������꣩
		conveyTool->SetUniform("useShadowTex", shaderData.bShadowTex, p);
		if (shaderData.bShadowTex)
		{
			conveyTool->SetTexture(shaderData.tShadowTex, 5, GL_TEXTURE5, "shadowTex", p);
		}
	}

	glDrawArrays(GL_TRIANGLES, 0, meshData.mesh.n_faces()*3);

	glBindVertexArray(0);
}
