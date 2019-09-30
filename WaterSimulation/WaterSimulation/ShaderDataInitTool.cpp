#include "ShaderDataInitTool.h"

ShaderDataInitTool* ShaderDataInitTool::instance = NULL;

void ShaderDataInitTool::InitVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexData, bool providedNormal, bool providedTex)
{
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
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
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);					//�Ȱ󶨣�����VAO��ֵʱ���ʹ��͵��ǵ�ǰ�󶨵�buffer

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

void ShaderDataInitTool::InitTextureWithFile(GLuint& texID, string texPath)
{
	glGenTextures(1, &texID);					//����һ������ID
	glBindTexture(GL_TEXTURE_2D, texID);		//��ʱ�󶨵���Ĭ������Ԫ0������֮��Ĵ����л�ָ���󶨵������ĸ���Ԫ
	//ָ����ͼ����
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//ͼƬ�ļ���ȡ
	int width, height;
	unsigned char* pResult = SOIL_load_image(texPath.c_str(), &width, &height, 0, SOIL_LOAD_RGB);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, pResult);
	//����һ��mipmap
	glGenerateMipmap(GL_TEXTURE_2D);
	//����󶨲��ͷ�
	glBindTexture(GL_TEXTURE_2D, 0);
	SOIL_free_image_data(pResult);
}
