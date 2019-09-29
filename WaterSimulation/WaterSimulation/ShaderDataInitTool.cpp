#include "ShaderDataInitTool.h"

ShaderDataInitTool* ShaderDataInitTool::instance = NULL;

void ShaderDataInitTool::SetVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexPos)
{
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	//创建顶点buffer
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);					//先绑定，在用VAO传值时，就传送的是当前绑定的buffer
	glBufferData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float), NULL, GL_STATIC_DRAW);
	//初始化数据
	glBufferSubData(GL_ARRAY_BUFFER, 0, vertexPos.size() * sizeof(float), &vertexPos[0]);
	//还需要用glVertexAttribPointer更新顶点属性指针
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	//解除绑定
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

void ShaderDataInitTool::SetVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexPos, vector<float>& vertexNormal)
{
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	//创建顶点buffer
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);					//先绑定，在用VAO传值时，就传送的是当前绑定的buffer
	glBufferData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float) + vertexNormal.size() * sizeof(float), NULL, GL_STATIC_DRAW);
	//初始化数据
	glBufferSubData(GL_ARRAY_BUFFER, 0, vertexPos.size() * sizeof(float), &vertexPos[0]);
	glBufferSubData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float), vertexNormal.size() * sizeof(float), &vertexNormal[0]);
	//还需要用glVertexAttribPointer更新顶点属性指针
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(vertexPos.size() * sizeof(float)));
	//解除绑定
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

void ShaderDataInitTool::SetVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexPos, vector<float>& vertexNormal, vector<float>& vertexTex)
{
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	//创建顶点buffer
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);					//先绑定，在用VAO传值时，就传送的是当前绑定的buffer
	glBufferData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float) + vertexNormal.size() * sizeof(float) + vertexTex.size() * sizeof(float), NULL, GL_STATIC_DRAW);
	//初始化数据
	glBufferSubData(GL_ARRAY_BUFFER, 0, vertexPos.size() * sizeof(float), &vertexPos[0]);
	glBufferSubData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float), vertexNormal.size() * sizeof(float), &vertexNormal[0]);
	glBufferSubData(GL_ARRAY_BUFFER, vertexPos.size() * sizeof(float) + vertexNormal.size() * sizeof(float), vertexTex.size() * sizeof(float), &vertexTex[0]);
	//还需要用glVertexAttribPointer更新顶点属性指针
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(vertexPos.size() * sizeof(float)));
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)((vertexPos.size() * sizeof(float)) + (vertexNormal.size() * sizeof(float))));
	//解除绑定
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

void ShaderDataInitTool::SetTexture(GLuint& texID, string texPath)
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
