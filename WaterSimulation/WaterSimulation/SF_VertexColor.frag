#version 430 

in vec3 posW;
in vec3 colorW;

out vec4 FragColor;

void main()
{
	FragColor = vec4(colorW, 1.0);
}