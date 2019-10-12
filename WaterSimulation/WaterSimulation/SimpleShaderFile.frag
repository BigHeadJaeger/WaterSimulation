#version 430 

in vec3 posW;
in vec3 normalW;
in vec2 TexCoord;

uniform vec3 color;

out vec4 FragColor;

void main()
{
	FragColor = vec4(color, 1.0);
}