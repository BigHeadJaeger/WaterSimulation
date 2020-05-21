#include "ShaderData.h"

// 控制物体整体的运动
void ShaderData::UpdateMatrix(Transform& t)
{
	world = translate(mat4(1.0), t.position);
	world = scale(world, t.scaler);
	if (t.rotation.x != 0)
		world = rotate(world, t.rotation.x, vec3(1.0, 0.0, 0.0));
	if (t.rotation.y != 0)
		world = rotate(world, t.rotation.y, vec3(0.0, 1.0, 0.0));
	if (t.rotation.z != 0)
		world = rotate(world, t.rotation.z, vec3(0.0, 0.0, 1.0));
	worldInvTranspose = transpose(inverse(world));
	worldViewProj = MainCamera::GetInstance()->pro * MainCamera::GetInstance()->view * world;
}

// VAO和VBO初始化、更新的接口
void ShaderData::InitVertexBuffer(vector<float>& vertexData, bool providedNormal, bool providedTex, GLenum usage)
{
	ShaderDataInitTool* tool = ShaderDataInitTool::GetShaderDataInitTool();
	tool->InitVertexBuffer(VAO, VBO, vertexData, providedNormal, providedTex, usage);
}
void ShaderData::UpdateVertexBuffer(vector<float>& vertexData, bool providedNormal, bool providedTex)
{
	ShaderDataInitTool* tool = ShaderDataInitTool::GetShaderDataInitTool();
	tool->UpdateVertexBuffer(VAO, VBO, vertexData, providedNormal, providedTex);
}

void UE4ShaderData::InitTexture(TEXTURETYPE type, string texPath)
{
	ShaderDataInitTool* tool = ShaderDataInitTool::GetShaderDataInitTool();
	switch (type)
	{
	case ALBEDO:
		tool->InitTextureWithFile(tAlbedo, texPath);
		break;
	case METALLIC:
		tool->InitTextureWithFile(tMetallic, texPath);
		break;
	case ROUGHNESS:
		tool->InitTextureWithFile(tRoughness, texPath);
		break;
	case AO:
		tool->InitTextureWithFile(tAo, texPath);
		break;
	case NORMAL:
		tool->InitTextureWithFile(tNormal, texPath);
		break;
	default:
		break;
	}
}

void VCShaderData::InitVertexBuffer(vector<float>& vertexData, GLenum usage)
{
	ShaderDataInitTool* tool = ShaderDataInitTool::GetShaderDataInitTool();
	tool->InitBufferVC(VAO, VBO, vertexData, usage);
}
