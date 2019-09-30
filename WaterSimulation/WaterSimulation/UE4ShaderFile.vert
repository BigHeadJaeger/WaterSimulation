#version 430

layout(location=0) in vec3 VertexPosition;
layout(location=1) in vec3 normal;
layout(location=2) in vec2 texCoord;


uniform mat4 world;
uniform mat4 worldViewProj;
uniform mat4 worldInvTranspose;

uniform mat4 depthBiasMVP;			//阴影贴图纹理修正矩阵


out vec3 posW;					//光照计算需要物体的世界坐标
out vec3 normalW;				//顶点法向量的世界坐标
out vec2 TexCoord;				//顶点纹理坐标

out vec4 shadowCoord;			//阴影贴图纹理坐标,同时z坐标也是当前点在光源视角下深度

void main()
{
	posW=(world*vec4(VertexPosition,1.0)).xyz;
	normalW=(worldInvTranspose*vec4(normal,1.0)).xyz;
	gl_Position=worldViewProj*vec4(VertexPosition,1.0);

	shadowCoord=depthBiasMVP*world*vec4(VertexPosition,1.0);	
	//TexCoord=shadowCoord.xy;
	TexCoord=texCoord;

}