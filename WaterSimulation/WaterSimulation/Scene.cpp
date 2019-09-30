#pragma once
#include "Scene.h"
#include<GLFW\glfw3.h>

void MyScene::Init()
{
	//��ʼ��glew
	glewInit();
	//��ʼ��Renderer�е�program
	UE4Renderer::GetRenderer()->InitProgram("test.vert", "test.frag");

	//pShadowTex.SetShader("shadowTex.v", "shadowTex.f");


	//ָ����Դ����
	lightPos = vec3(-3, 5, 3);
	lightColor = vec3(400.f, 400.f, 400.f);

	//��ʼ�������
	mainCamera = new Camera();
	mainCamera->Init(vec3(0, 2, 4), vec3(0, 0, 0));

	//SetDrawMode(drawMode.isLine, false);
	drawMode.isLine = false;

	
	

	//ָ������PBR����
	MeshObject* cow = new MeshObject();
	cow->SetName("cow");
	cow->readObjFile("OBJ\\cow.obj");
	cow->SetRenderer(UE4RENDERER);
	cow->InitBufferData();
	cow->GetTransform().SetPosition(vec3(0, 0, 0));
	cow->GetTransform().SetScaler(vec3(3.0));
	dynamic_cast<UE4ShaderData*>(cow->GetShaderData())->InitTexture(ALBEDO, "");

	objects.insert(pair<string, Object*>(cow->GetName(), cow));


	//myBox.InitDirectBox(1, 1, 1);					//���㡢������Ϣ��ʼ��
	//myBox.InitBuffers();							//�����ʼ��
	//myBox.CoorDataInit(vec3(-0.5, 0.5, 1.0)/*, vec3(3.0, 0.1, 3.0)*/);
	//myBox.TextureUseInfo(true, true, true, true, true, true);				//����ʹ����Ϣ��ʼ��
	//myBox.InitTexture(myBox.TAlbedo, "Material\\metalgrid2-dx\\metalgrid2_basecolor.png");		//������Ҫ������
	//myBox.InitTexture(myBox.TNormal, "Material\\metalgrid2-dx\\metalgrid2_normal-dx.png");
	//myBox.InitTexture(myBox.TAo, "Material\\metalgrid2-dx\\metalgrid2_AO.png");
	//myBox.InitTexture(myBox.TRoughness, "Material\\metalgrid2-dx\\metalgrid2_roughness.png");
	//myBox.InitTexture(myBox.TMetallic, "Material\\metalgrid2-dx\\metalgrid2_metallic.png");
	//myBox.LightUseInfo(true);						//�����Ƿ���չ���
	//myBox.InitMaterial(vec4(0.5f, 0.5f, 0.5f, 1.0f), vec4(0.5, 0.5, 0.5, 1.0f), vec4(0.7f, 0.7f, 0.7f, 1.0f), 7.0f);


	//STInit();
}

void MyScene::InitKeys()
{
	keys.insert(pair<KEYNAME, Key>(BTNW, Key(BTNW)));
	keys.insert(pair<KEYNAME, Key>(BTNA, Key(BTNA)));
	keys.insert(pair<KEYNAME, Key>(BTNS, Key(BTNS)));
	keys.insert(pair<KEYNAME, Key>(BTND, Key(BTND)));
	keys.insert(pair<KEYNAME, Key>(BTN1, Key(BTN1)));
	//keys.push_back(Key(BTNW));
}

void MyScene::Update()
{

	//�����ӽǾ���
	mainCamera->SetView();
	//����ͶӰ����
	mainCamera->SetPro();

	//��������object���¾���
	map<string, Object*>::iterator objs_it;
	for (objs_it = objects.begin(); objs_it != objects.end(); objs_it++)
	{
		(*objs_it).second->Update();
	}

	//cow.UpdateMatrix(mainCamera);
	//myBox.SetObjMat(camera.view, camera.pro);
	//myBucket.SetObjMat(camera.view, camera.pro);
	//myGrid.SetObjMat(camera.view, camera.pro);

	//��������key����ִ��key��ǰ�󶨵��¼�
	map<KEYNAME, Key>::iterator keys_it;
	for (keys_it = keys.begin(); keys_it != keys.end(); keys_it++)
	{
		keys_it->second.Execute();
	}
}

void MyScene::Draw()
{

	//NormalShadowMap();

	//����ÿһ������ʱ���Ƚ����еĸ������������ؾ��󣬲��ʵ���Ϣ����shader��


	//����,��������������գ�����Ч���Ŀ�����blend�� cull֮��ģ�
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);					//��ά����Ҫ���������޳�

	if (drawMode.isLine)
	{
		glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else
	{
		glEnable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	//glFrontFace(GL_CW);
	//glCullFace(GL_FRONT);

	//glUseProgram(p1.p);						//������ɫ������

	map<string, Object*>::iterator objs_it;
	for (objs_it = objects.begin(); objs_it != objects.end(); objs_it++)
	{
		(*objs_it).second->Draw();
	}

}