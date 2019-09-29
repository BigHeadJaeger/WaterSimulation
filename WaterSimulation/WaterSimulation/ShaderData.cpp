#include "ShaderData.h"

void ShaderData::UpdateMatrix(Transform& t, Camera& camera)
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
	worldViewProj = camera.pro * camera.view * world;
}

void ShaderData::InitTexture(GLuint& texID, string texPath)
{
	glGenTextures(1, &texID);					//生成一个纹理ID
	glBindTexture(GL_TEXTURE_2D, texID);		//此时绑定到了默认纹理单元0处，在之后的代码中会指定绑定到具体哪个单元
	//指定贴图方法
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//图片文件读取
	int width, height;
	unsigned char* pResult = SOIL_load_image(texPath.c_str(), &width, &height, 0, SOIL_LOAD_RGB);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, pResult);
	//生成一个mipmap
	glGenerateMipmap(GL_TEXTURE_2D);
	//解除绑定并释放
	glBindTexture(GL_TEXTURE_2D, 0);
	SOIL_free_image_data(pResult);
}

void ShaderData::InitVertexBuffer(MeshData& meshData)
{
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);

	//生成顶点buffer(VBO)并绑定
	//MyGeometryGenerator::GetInstance()->CreateBox(20.0f, 0.3f, 20.0f, box);
	//转换到一维数组中
	std::vector<float> vertexData;
	vector<float> vertexPos;
	vector<float> vertexNormal;
	vector<float> vertexTex;
	for (Mesh::FaceIter f_it = meshData.mesh.faces_begin(); f_it != meshData.mesh.faces_end(); f_it++)
	{
		for (Mesh::FaceVertexIter fv_ccwit = meshData.mesh.fv_iter(*f_it); fv_ccwit.is_valid(); fv_ccwit++)
		{
			vertexPos.push_back(meshData.mesh.point(*fv_ccwit).data()[0]);
			vertexPos.push_back(meshData.mesh.point(*fv_ccwit).data()[1]);
			vertexPos.push_back(meshData.mesh.point(*fv_ccwit).data()[2]);

			vertexNormal.push_back(meshData.mesh.normal(*fv_ccwit).data()[0]);
			vertexNormal.push_back(meshData.mesh.normal(*fv_ccwit).data()[1]);
			vertexNormal.push_back(meshData.mesh.normal(*fv_ccwit).data()[2]);

			if (meshData.providedTex)
			{
				vertexTex.push_back(meshData.mesh.texcoord2D(*fv_ccwit).data()[0]);
				vertexTex.push_back(meshData.mesh.texcoord2D(*fv_ccwit).data()[1]);
			}
		}
	}

	//创建顶点buffer
	glGenBuffers(1, &VertexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBuffer);					//先绑定，在用VAO传值时，就传送的是当前绑定的buffer
	//开辟空间
	if (meshData.providedTex)
		glBufferData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float) + vertexNormal.size() * sizeof(float) + vertexTex.size() * sizeof(float), NULL, GL_STATIC_DRAW);
	else
		glBufferData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float) + vertexNormal.size() * sizeof(float), NULL, GL_STATIC_DRAW);

	//初始化数据
	glBufferSubData(GL_ARRAY_BUFFER, 0, vertexPos.size() * sizeof(float), &vertexPos[0]);
	glBufferSubData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float), vertexNormal.size() * sizeof(float), &vertexNormal[0]);
	if (meshData.providedTex)
	{
		glBufferSubData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float) + vertexNormal.size() * sizeof(float), vertexTex.size() * sizeof(float), &vertexTex[0]);
	}

	//还需要用glVertexAttribPointer更新顶点属性指针
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(vertexPos.size() * sizeof(float)));
	if (meshData.providedTex)
	{
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)((vertexPos.size() * sizeof(float)) + (vertexNormal.size() * sizeof(float))));
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

}
