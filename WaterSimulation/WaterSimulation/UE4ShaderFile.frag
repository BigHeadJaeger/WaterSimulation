#version 430 

//这里采用的直接光照，不考虑环境光照，因此每一个点的入射光只考虑存在的光源个数，将总辐照度相加即可，如果考虑全局光照，则需要解积分计算

in vec3 posW;
in vec3 normalW;
in vec2 TexCoord;
in vec4 shadowCoord;


uniform bool useTexture;

//阴影贴图纹理采样器
uniform sampler2DShadow shadowTex;
uniform bool useShadowTex;

//物体PBR材质
//当前属性的值是否用纹理来决定
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

//光照信息(此处使用的是点光源)
uniform vec3 lightPos;
uniform vec3 lightColor; 

//眼睛位置
uniform vec3 eyePos;

const float PI=3.14159265359;
 
out vec4 FragColor;

//正态分布函数
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

//几何函数
float GeometrySchlickGGX(vec3 N,vec3 v,float roughness)
{
	//计算基于直接光照的k值
	float k=(roughness+1)*(roughness+1)/8.0;

	float Ndotv=max(dot(N,v),0.0);
	float nom=Ndotv;

	float denom=Ndotv*(1.0-k)+k;

	return nom/denom;


}
float GeometrySmith(vec3 N,vec3 L,vec3 E,float roughness)			//几何阴影和几何遮蔽
{
	float GSGGXL=GeometrySchlickGGX(N,L,roughness);
	float GSGGXE=GeometrySchlickGGX(N,E,roughness);

	return GSGGXL*GSGGXE;
}

//菲涅尔方程
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

		//从各种贴图中获取数据
		if(useAlbedo)
			albedo=pow(texture(albedoMap,TexCoord).rgb,vec3(2.2));			//反射率纹理一般创建在rgb空间，所以需要转换到线性空间
		 
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
	//采用法线贴图
	if(useNormal)
		N = GetNormalFromMap();

    vec3 V = normalize(eyePos - posW);
	 
	vec3 F0=vec3(0.04);
	F0=mix(F0,albedo,metallic);

	vec3 reflectResult=vec3(0.0);


	//计算辐射度(根据入射方向以及夹角求)
	vec3 lightDir=normalize(lightPos-posW);				//计算光照向量
	vec3 H=normalize(V+lightDir);									//计算中间向量


	float distance=length(lightPos-posW);
	float attenuation=1.0/(distance*distance);			//计算衰减
	vec3 radiance=lightColor*attenuation;

	//计算BDRF中的镜面反射
	float D=DistributionGGX(H,N,roughness);
	float G=GeometrySmith(N,lightDir,V,roughness);
	vec3 F=Fresnel(clamp(dot(H,V),0.0,1.0),F0);

	vec3 nom=D*G*F;
	float denom=4 * max(dot(N, V), 0.0) * max(dot(N, lightDir), 0.0);	//0.001防止除0
	vec3 specular=nom/max(denom,0.001);

	//计算BDRF中的漫反射
	vec3 KS=F;					//菲涅尔中已经求出了反射光线的占比
	vec3 KD=vec3(1.0)-KS;
	KD*=1.0-metallic;			//因为金属没有漫反射，所以根据金属度重新决定漫反射


	float NdotL=max(dot(N,lightDir),0.0);
	reflectResult+=(KD*albedo/PI+specular)*radiance*NdotL;


	vec3 ambient=vec3(0.03)*albedo*ao;
	color=ambient+reflectResult;

	//伽马校正
	color = color / (color + vec3(1.0));
    color = pow(color, vec3(1.0/2.2));


	//如果当前点的深度大于此点在阴影贴图中的深度，说明这个点在阴影中
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