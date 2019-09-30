#version 430 

//������õ�ֱ�ӹ��գ������ǻ������գ����ÿһ����������ֻ���Ǵ��ڵĹ�Դ���������ܷ��ն���Ӽ��ɣ��������ȫ�ֹ��գ�����Ҫ����ּ���

in vec3 posW;
in vec3 normalW;
in vec2 TexCoord;
in vec4 shadowCoord;


uniform bool useTexture;

//��Ӱ��ͼ���������
uniform sampler2DShadow shadowTex;
uniform bool useShadowTex;

//����PBR����
//��ǰ���Ե�ֵ�Ƿ�������������
uniform sampler2D albedoMap;
uniform bool useAlbedo;
uniform vec3 albedoN;

uniform sampler2D normalMap;
uniform bool useNormal;
uniform vec3 normalN;

uniform sampler2D metallicMap;
uniform bool useMetallic;
uniform float metallicN;

uniform sampler2D roughnessMap;
uniform bool useRoughness;
uniform float roughnessN;

uniform sampler2D aoMap;
uniform bool useAO;
uniform float aoN;

//������Ϣ(�˴�ʹ�õ��ǵ��Դ)
uniform vec3 lightPos;
uniform vec3 lightColor; 

//�۾�λ��
uniform vec3 eyePos;

const float PI=3.14159265359;
 
out vec4 FragColor;

//��̬�ֲ�����
float DistributionGGX(vec3 H,vec3 N,float roughness)
{
	float a=roughness*roughness;
	float a2=a*a;
	float NdotH=max(dot(N,H),0.0);
	float NdotH2=NdotH*NdotH;

	float denom=PI*(NdotH2*(a2-1)+1)*(NdotH2*(a2-1)+1);
	float nom=a2;

	return nom/max(denom,0.001);
}

//���κ���
float GeometrySchlickGGX(vec3 N,vec3 v,float roughness)
{
	//�������ֱ�ӹ��յ�kֵ
	float k=(roughness+1)*(roughness+1)/8.0;

	float Ndotv=max(dot(N,v),0.0);
	float nom=Ndotv;

	float denom=Ndotv*(1.0-k)+k;

	return nom/denom;


}
float GeometrySmith(vec3 N,vec3 L,vec3 E,float roughness)			//������Ӱ�ͼ����ڱ�
{
	float GSGGXL=GeometrySchlickGGX(N,L,roughness);
	float GSGGXE=GeometrySchlickGGX(N,E,roughness);

	return GSGGXL*GSGGXE;
}

//����������
vec3 Fresnel(float cosTheta,vec3 F0)
{
	return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 GetNormalFromMap()
{
	vec3 tangentNormal = texture(normalMap, TexCoord).xyz * 2.0 - 1.0;

    vec3 Q1  = dFdx(posW);
    vec3 Q2  = dFdy(posW);
    vec2 st1 = dFdx(TexCoord);
    vec2 st2 = dFdy(TexCoord);

    vec3 N   = normalize(normalW);
    vec3 T  = normalize(Q1*st2.t - Q2*st1.t);
    vec3 B  = -normalize(cross(N, T));
    mat3 TBN = mat3(T, B, N);

    return normalize(TBN * tangentNormal);
}

void main() 
{
	vec3 color=vec3(1.0,0.0,0.0);
	if(useTexture)
	{
		vec3 albedo;
		float roughness;
		float ao;
		float metallic;

		//�Ӹ�����ͼ�л�ȡ����
		if(useAlbedo)
			albedo=pow(texture(albedoMap,TexCoord).rgb,vec3(2.2));			//����������һ�㴴����rgb�ռ䣬������Ҫת�������Կռ�
		 
	if(useMetallic)
		metallic=texture(metallicMap,TexCoord).r;
	else
		metallic=metallicN;


	if(useRoughness)
		roughness=texture(roughnessMap,TexCoord).r;
	else
		roughness=roughnessN;

	if(useAO)
		ao=texture(aoMap,TexCoord).r;
	else
		ao=aoN;

		

	vec3 N;
	//���÷�����ͼ
	if(useNormal)
		N = GetNormalFromMap();

    vec3 V = normalize(eyePos - posW);
	 
	vec3 F0=vec3(0.04);
	F0=mix(F0,albedo,metallic);

	vec3 reflectResult=vec3(0.0);


	//��������(�������䷽���Լ��н���)
	vec3 lightDir=normalize(lightPos-posW);				//�����������
	vec3 H=normalize(V+lightDir);									//�����м�����


	float distance=length(lightPos-posW);
	float attenuation=1.0/(distance*distance);			//����˥��
	vec3 radiance=lightColor*attenuation;

	//����BDRF�еľ��淴��
	float D=DistributionGGX(H,N,roughness);
	float G=GeometrySmith(N,lightDir,V,roughness);
	vec3 F=Fresnel(clamp(dot(H,V),0.0,1.0),F0);

	vec3 nom=D*G*F;
	float denom=4 * max(dot(N, V), 0.0) * max(dot(N, lightDir), 0.0);	//0.001��ֹ��0
	vec3 specular=nom/max(denom,0.001);

	//����BDRF�е�������
	vec3 KS=F;					//���������Ѿ�����˷�����ߵ�ռ��
	vec3 KD=vec3(1.0)-KS;
	KD*=1.0-metallic;			//��Ϊ����û�������䣬���Ը��ݽ��������¾���������


	float NdotL=max(dot(N,lightDir),0.0);
	reflectResult+=(KD*albedo/PI+specular)*radiance*NdotL;


	vec3 ambient=vec3(0.03)*albedo*ao;
	color=ambient+reflectResult;

	//٤��У��
	color = color / (color + vec3(1.0));
    color = pow(color, vec3(1.0/2.2));


	//�����ǰ�����ȴ��ڴ˵�����Ӱ��ͼ�е���ȣ�˵�����������Ӱ��
	float bias=0.000;
	float visibility=1.0;
	if(useShadowTex)
	{
		//if(texture(shadowTex,shadowCoord.xy).z+bias<shadowCoord.z)
		if(texture(shadowTex,vec3(shadowCoord.xy,shadowCoord.z-bias))!=1)
		visibility=0.5;	
	}
	color*=visibility;
	}



    FragColor = vec4(color, 1.0); 
	//FragColor = vec4(1.0,0.0,0.0, 1.0); 
} 