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
	vector<float> vertexPos;
	vector<float> vertexNormal;
	vector<float> vertexTex;
	for (int i = 0; i < vertexData.size(); i = i + 8)
	{
		vertexPos.push_back(vertexData[i]);
		vertexPos.push_back(vertexData[i + 1]);
		vertexPos.push_back(vertexData[i + 2]);
		if (providedNormal)
		{
			vertexPos.push_back(vertexData[i + 3]);
			vertexPos.push_back(vertexData[i + 4]);
			vertexPos.push_back(vertexData[i + 5]);
		}
		if (providedTex)
		{
			vertexPos.push_back(vertexData[i + 6]);
			vertexPos.push_back(vertexData[i + 7]);
		}
	}

	//��������buffer
	glGenBuffers(1, &VertexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBuffer);					//�Ȱ󶨣�����VAO��ֵʱ���ʹ��͵��ǵ�ǰ�󶨵�buffer
	
	GLsizeiptr sumSize = 0;
	sumSize += vertexPos.size() * sizeof(float);
	if (providedNormal)
		sumSize += vertexNormal.size() * sizeof(float);
	if (providedTex)
		sumSize += vertexTex.size() * sizeof(float);
	//���ٿռ�
	glBufferData(GL_ARRAY_BUFFER, sumSize, NULL, GL_STATIC_DRAW);

	GLintptr offset = 0;
	//��ʼ������
	glBufferSubData(GL_ARRAY_BUFFER, offset, vertexPos.size() * sizeof(float), &vertexPos[0]);
	offset += vertexPos.size() * sizeof(float);
	if (providedNormal)
	{
		glBufferSubData(GL_ARRAY_BUFFER, offset, vertexNormal.size() * sizeof(float), &vertexNormal[0]);
		offset += vertexNormal.size() * sizeof(float);
	}
	if (providedTex)
	{
		glBufferSubData(GL_ARRAY_BUFFER, offset, vertexTex.size() * sizeof(float), &vertexTex[0]);
		offset += vertexNormal.size() * sizeof(float);
	}

	GLint index = 0;
	GLsizeiptr size = 0;
	//����Ҫ��glVertexAttribPointer���¶�������ָ��
	glEnableVertexAttribArray(index);
	glVertexAttribPointer(index, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)size);
	index++;
	size += vertexPos.size() * sizeof(float);
	if (providedNormal)
	{
		glEnableVertexAttribArray(index);
		glVertexAttribPointer(index, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)size);
		index++;
		size += vertexNormal.size() * sizeof(float);
	}
	if (providedTex)
	{
		glEnableVertexAttribArray(index);
		glVertexAttribPointer(index, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)size);
		index++;
		size += vertexNormal.size() * sizeof(float);
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}
