#version 430

layout(location=0) in vec3 VertexPosition;
layout(location=1) in vec3 color;



uniform mat4 world;
uniform mat4 worldViewProj;
uniform mat4 worldInvTranspose;

out vec3 posW;					//片段着色器物体世界坐标
out vec3 colorW;					//片段着色器物体顶点颜色

void main()
{
	posW=(world*vec4(VertexPosition,1.0)).xyz;
	colorW = color;
	//normalW=(worldInvTranspose*vec4(normal,1.0)).xyz;
	gl_Position=worldViewProj*vec4(VertexPosition,1.0);

	//shadowCoord=depthBiasMVP*world*vec4(VertexPosition,1.0);	
	//TexCoord=shadowCoord.xy;
	//TexCoord=texCoord;

}