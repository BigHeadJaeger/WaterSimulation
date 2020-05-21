#include "ShaderDataInitTool.h"

ShaderDataInitTool* ShaderDataInitTool::instance = NULL;



void ShaderDataInitTool::InitVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexData, bool providedNormal, bool providedTex, GLenum usage)
{
	glDeleteVertexArrays(1, &VAO);
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	vector<float> vertexPos;
	vector<float> vertexNormal;
	vector<float> vertexTex;
	VertexDataUnpack(vertexData, vertexPos, vertexNormal, vertexTex);

	//创建顶点buffer
	glDeleteBuffers(1, &VBO);
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);					//先绑定，在用VAO传值时，就传送的是当前绑定的buffer

	GLsizeiptr sumSize = 0;
	sumSize += vertexPos.size() * sizeof(float);
	if (providedNormal)
		sumSize += vertexNormal.size() * sizeof(float);
	if (providedTex)
		sumSize += vertexTex.size() * sizeof(float);
	//开辟空间
	glBufferData(GL_ARRAY_BUFFER, sumSize, NULL, usage);

	GLintptr offset = 0;
	//初始化数据
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
	//还需要用glVertexAttribPointer更新顶点属性指针
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

void ShaderDataInitTool::UpdateVertexBuffer(GLuint& VAO, GLuint& VBO, vector<float>& vertexData, bool providedNormal, bool providedTex)
{
	if (!glIsVertexArray(VAO))
	{
		cout << "VAO not exist, need init before update" << endl;
		return;
	}
	if (!glIsBuffer(VBO))
	{
		cout << "VBO not exist, need init before update" << endl;
		return;
	}

	glBindVertexArray(VAO);
	vector<float> vertexPos;
	vector<float> vertexNormal;
	vector<float> vertexTex;
	VertexDataUnpack(vertexData, vertexPos, vertexNormal, vertexTex);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	GLintptr offset = 0;
	// 重新分配新的数据
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
	//还需要用glVertexAttribPointer更新顶点属性指针
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

void ShaderDataInitTool::InitBufferVC(GLuint& VAO, GLuint& VBO, vector<float>& vertexData, GLenum usage)
{
	glDeleteVertexArrays(1, &VAO);
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	vector<float> vertexPos;
	vector<float> vertexColor;
	for (int i = 0; i < vertexData.size(); i = i + 6)
	{
		vertexPos.push_back(vertexData[i]);
		vertexPos.push_back(vertexData[i + 1]);
		vertexPos.push_back(vertexData[i + 2]);
		vertexColor.push_back(vertexData[i + 3]);
		vertexColor.push_back(vertexData[i + 4]);
		vertexColor.push_back(vertexData[i + 5]);
	}

	//创建顶点buffer
	glDeleteBuffers(1, &VBO);
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);					//先绑定，在用VAO传值时，就传送的是当前绑定的buffer

	GLsizeiptr sumSize = 0;
	sumSize += vertexPos.size() * sizeof(float);
	sumSize += vertexColor.size() * sizeof(float);

	//开辟空间
	glBufferData(GL_ARRAY_BUFFER, sumSize, NULL, usage);

	GLintptr offset = 0;
	//初始化数据
	glBufferSubData(GL_ARRAY_BUFFER, offset, vertexPos.size() * sizeof(float), &vertexPos[0]);
	offset += vertexPos.size() * sizeof(float);

	glBufferSubData(GL_ARRAY_BUFFER, offset, vertexColor.size() * sizeof(float), &vertexColor[0]);
	offset += vertexColor.size() * sizeof(float);

	GLint index = 0;
	GLsizeiptr size = 0;
	//还需要用glVertexAttribPointer更新顶点属性指针
	glEnableVertexAttribArray(index);
	glVertexAttribPointer(index, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)size);
	index++;
	size += vertexPos.size() * sizeof(float);

	glEnableVertexAttribArray(index);
	glVertexAttribPointer(index, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)size);
	index++;
	size += vertexColor.size() * sizeof(float);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}


void ShaderDataInitTool::VertexDataUnpack(vector<float>& vertexData, vector<float>& vertexPos, vector<float>& vertexNormal, vector<float>& vertexTex)
{
	for (int i = 0; i < vertexData.size(); i = i + 8)
	{
		vertexPos.push_back(vertexData[i]);
		vertexPos.push_back(vertexData[i + 1]);
		vertexPos.push_back(vertexData[i + 2]);
		vertexNormal.push_back(vertexData[i + 3]);
		vertexNormal.push_back(vertexData[i + 4]);
		vertexNormal.push_back(vertexData[i + 5]);
		vertexTex.push_back(vertexData[i + 6]);
		vertexTex.push_back(vertexData[i + 7]);
	}
}
