#include "Program.h"

void ShaderProgram::printShaderInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten = 0;
	char *infoLog;

	glGetShaderiv(obj, GL_INFO_LOG_LENGTH, &infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
		printf("%s\n", infoLog);
		free(infoLog);
	}
}

void ShaderProgram::printProgramInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten = 0;
	char *infoLog;

	glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
		printf("%s\n", infoLog);
		free(infoLog);
	}
}

char * ShaderProgram::textFileRead(const char * fn)
{  //read the shader code
	FILE *fp;
	char *content = NULL;
	int count = 0;

	if (fn != NULL) {
		fopen_s(&fp, fn, "rt");

		if (fp != NULL) {
			fseek(fp, 0, SEEK_END);
			count = ftell(fp);
			rewind(fp);

			if (count > 0) {
				content = (char *)malloc(sizeof(char) * (count + 1));
				count = fread(content, sizeof(char), count, fp);
				content[count] = '\0';
			}
			fclose(fp);
		}
		else
			std::cout << "Fail to open the shader file" << std::endl;
	}
	return content;
}

void ShaderProgram::SetShader(const char*VSFile, const char*FSFile)
{
	char*vs, *fs;
	//���ݲ���ȷ������ʲô��ɫ��
	v = glCreateShader(GL_VERTEX_SHADER);
	f = glCreateShader(GL_FRAGMENT_SHADER);
	//��ȡshader���ַ���
	vs = textFileRead(VSFile);
	fs = textFileRead(FSFile);

	const char*vv = vs;
	const char*ff = fs;
	//Ϊ������bufferָ��shaderԴ
	glShaderSource(v, 1, &vv, NULL);
	glShaderSource(f, 1, &ff, NULL);

	free(vs);
	free(fs);
	//����shader
	glCompileShader(v);
	glCompileShader(f);

	printShaderInfoLog(v);
	printShaderInfoLog(f);
	//����һ��shader����
	p = glCreateProgram();
	//��������ɫ���󶨵�shader������
	glAttachShader(p, v);
	glAttachShader(p, f);
	//����shader����OpenGL������
	glLinkProgram(p);


	printProgramInfoLog(p);

}
