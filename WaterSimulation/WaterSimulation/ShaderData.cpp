#include "ShaderData.h"

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
	worldViewProj = mainCamera->pro * mainCamera->view * world;
}

void ShaderData::InitVertexBuffer(vector<float>& vertexData, bool providedNormal, bool providedTex)
{
	ShaderDataInitTool* tool = ShaderDataInitTool::GetShaderDataInitTool();
	tool->InitVertexBuffer(VAO, VBO, vertexData, providedNormal, providedTex);
}
