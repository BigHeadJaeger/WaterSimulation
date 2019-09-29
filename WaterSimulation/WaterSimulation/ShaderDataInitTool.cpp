#include "ShaderDataInitTool.h"

ShaderDataInitTool* ShaderDataInitTool::instance = NULL;

void ShaderDataInitTool::SetVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexPos)
{
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	//��������buffer
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);					//�Ȱ󶨣�����VAO��ֵʱ���ʹ��͵��ǵ�ǰ�󶨵�buffer
	glBufferData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float), NULL, GL_STATIC_DRAW);
	//��ʼ������
	glBufferSubData(GL_ARRAY_BUFFER, 0, vertexPos.size() * sizeof(float), &vertexPos[0]);
	//����Ҫ��glVertexAttribPointer���¶�������ָ��
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	//�����
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

void ShaderDataInitTool::SetVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexPos, vector<float>& vertexNormal)
{
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	//��������buffer
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);					//�Ȱ󶨣�����VAO��ֵʱ���ʹ��͵��ǵ�ǰ�󶨵�buffer
	glBufferData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float) + vertexNormal.size() * sizeof(float), NULL, GL_STATIC_DRAW);
	//��ʼ������
	glBufferSubData(GL_ARRAY_BUFFER, 0, vertexPos.size() * sizeof(float), &vertexPos[0]);
	glBufferSubData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float), vertexNormal.size() * sizeof(float), &vertexNormal[0]);
	//����Ҫ��glVertexAttribPointer���¶�������ָ��
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(vertexPos.size() * sizeof(float)));
	//�����
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

void ShaderDataInitTool::SetVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexPos, vector<float>& vertexNormal, vector<float>& vertexTex)
{
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	//��������buffer
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);					//�Ȱ󶨣�����VAO��ֵʱ���ʹ��͵��ǵ�ǰ�󶨵�buffer
	glBufferData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float) + vertexNormal.size() * sizeof(float) + vertexTex.size() * sizeof(float), NULL, GL_STATIC_DRAW);
	//��ʼ������
	glBufferSubData(GL_ARRAY_BUFFER, 0, vertexPos.size() * sizeof(float), &vertexPos[0]);
	glBufferSubData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float), vertexNormal.size() * sizeof(float), &vertexNormal[0]);
	glBufferSubData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float) + vertexNormal.size() * sizeof(float), vertexTex.size() * sizeof(float), &vertexTex[0]);
	//����Ҫ��glVertexAttribPointer���¶�������ָ��
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(vertexPos.size() * sizeof(float)));
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)((vertexPos.size() * sizeof(float)) + (vertexNormal.size() * sizeof(float))));
	//�����
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

void ShaderDataInitTool::SetTexture(GLuint& texID, string texPath)
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
