#version 430

layout(location=0) in vec3 VertexPosition;
layout(location=1) in vec3 normal;
layout(location=2) in vec2 texCoord;


uniform mat4 world;
uniform mat4 worldViewProj;
uniform mat4 worldInvTranspose;

uniform mat4 depthBiasMVP;			//��Ӱ��ͼ������������


out vec3 posW;					//���ռ�����Ҫ�������������
out vec3 normalW;				//���㷨��������������
out vec2 TexCoord;				//������������

out vec4 shadowCoord;			//��Ӱ��ͼ��������,ͬʱz����Ҳ�ǵ�ǰ���ڹ�Դ�ӽ������

void main()
{
	posW=(world*vec4(VertexPosition,1.0)).xyz;
	normalW=(worldInvTranspose*vec4(normal,1.0)).xyz;
	gl_Position=worldViewProj*vec4(VertexPosition,1.0);

	shadowCoord=depthBiasMVP*world*vec4(VertexPosition,1.0);	
	//TexCoord=shadowCoord.xy;
	TexCoord=texCoord;

}