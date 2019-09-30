#pragma once
#define GLEW_STATIC
#include<GL/glew.h>
#include<iostream>
class ShaderProgram
{
public:
	GLuint p;
	GLuint v;
	GLuint f;
public:
	void SetShader(const char*VSFile, const char*FSFile);
private:
	void printShaderInfoLog(GLuint obj);			//shader编译错误信息输出
	void printProgramInfoLog(GLuint obj);			//shader程序链接错误信息输出
	char *textFileRead(const char *fn);
};