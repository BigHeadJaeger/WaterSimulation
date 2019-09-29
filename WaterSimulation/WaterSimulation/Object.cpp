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

	glBindVertexArray(shaderData.VAO);				//绑定前面设置好的VAO
									//传递坐标变换矩阵
	conveyTool->SetUniform("worldViewProj", shaderData.worldViewProj, p);
	conveyTool->SetUniform("world", shaderData.world, p);
	conveyTool->SetUniform("worldInvTranspose", shaderData.worldInvTranspose, p);

	//SetUniform("depthBiasMVP", depthBiasMVP, p);


	//根据参数上对纹理的选择，将需要的纹理传入着色器
	//先将是否使用纹理传入shader
	conveyTool->SetUniform("useTexture", shaderData.bUseTexture, p);
	if (shaderData.bUseTexture)
	{
		//当物体使用纹理且有纹理坐标
		if (meshData.providedTex)
		{
			//基础反射贴图
			conveyTool->SetUniform("useAlbedo", shaderData.bAlbedo, p);
			if (shaderData.bAlbedo)
			{
				conveyTool->SetTexture(shaderData.tAlbedo, 0, GL_TEXTURE0, "albedoMap", p);

			}

			//法线贴图
			conveyTool->SetUniform("useNormal", shaderData.bNormal, p);
			if (shaderData.bNormal)
			{
				conveyTool->SetTexture(shaderData.tNormal, 1, GL_TEXTURE1, "normalMap", p);

			}

			//金属度贴图
			conveyTool->SetUniform("useMetallic", shaderData.bMetallic, p);
			if (shaderData.bMetallic)
			{
				conveyTool->SetTexture(shaderData.tMetallic, 2, GL_TEXTURE2, "metallicMap", p);

			}
			else
			{
				//此处暂时直接将没有金属贴图的金属度用数字传入shader
				conveyTool->SetUniform("metallicN", 0.2f, p);
			}

			//粗糙贴图
			conveyTool->SetUniform("useRoughness", shaderData.bRoughness, p);
			if (shaderData.bRoughness)
			{
				conveyTool->SetTexture(shaderData.tRoughness, 3, GL_TEXTURE3, "roughnessMap", p);

			}

			//环境光ao贴图
			conveyTool->SetUniform("useAO", shaderData.bAo, p);
			if (shaderData.bAo)
			{
				conveyTool->SetTexture(shaderData.tAo, 4, GL_TEXTURE4, "aoMap", p);

			}


		}

		//阴影贴图（与其它的纹理不是同一种类型，不需要物体自身的纹理坐标）
		conveyTool->SetUniform("useShadowTex", shaderData.bShadowTex, p);
		if (shaderData.bShadowTex)
		{
			conveyTool->SetTexture(shaderData.tShadowTex, 5, GL_TEXTURE5, "shadowTex", p);
		}
	}

	glDrawArrays(GL_TRIANGLES, 0, meshData.mesh.n_faces()*3);

	glBindVertexArray(0);
}
