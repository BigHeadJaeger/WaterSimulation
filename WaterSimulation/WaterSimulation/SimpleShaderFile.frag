in vec3 posW;
in vec3 normalW;
in vec2 TexCoord;

uniform vec3 color;

void main()
{
	FragColor = vec4(color, 1.0);
}