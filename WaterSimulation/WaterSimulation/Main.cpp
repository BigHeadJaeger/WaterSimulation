#pragma once
#include"Scene.h"
#include"Const.h"
#include<GLFW\glfw3.h>

//const GLuint WIDTH = 1200, HEIGHT = 1000;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);

MyScene scene;

float deltaTime = 0.0f;
float lastFrame = 0.0f;

int main(void)
{
	//��������
	GLFWwindow* window;

	//��ʼ��glfw
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(WIDTH, HEIGHT, "Hello World", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* �ú���ָ���ڵ�ǰ�����߳��ϵĴ�����Ⱦ����ΪOpenGL��OpenGL ES���� */
	glfwMakeContextCurrent(window);

	//ָ���¼�
	glfwSetKeyCallback(window, key_callback);			//�����¼�
	//glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	glfwSetMouseButtonCallback(window, mouse_button_callback);		//��갴���¼�

	glfwSetCursorPosCallback(window, cursor_position_callback);		//���ָ���¼�


	//��������

	scene.Init();


	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		//����һ֡�ļ��
		float currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;



		//�������ºͻ���
		scene.Update(deltaTime);
		scene.Draw();

		//����˫����
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
		//glfwWaitEvents();
	}

	glfwTerminate();
	return 0;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	//��keep�¼���һ��Ҫ��up�¼�
	if (key == GLFW_KEY_W && action == GLFW_PRESS)
		scene.keys[BTNW].BindKeepEvent([]() {
		MainCamera::GetInstance()->Walk(MainCamera::GetInstance()->cameraSpeed * deltaTime);
			});
	if (key == GLFW_KEY_W && action == GLFW_RELEASE)
		scene.keys[BTNW].BindUpEvent([]() {
			});

	if (key == GLFW_KEY_S && action == GLFW_PRESS)
		scene.keys[BTNS].BindKeepEvent([]() {
		MainCamera::GetInstance()->Walk(-MainCamera::GetInstance()->cameraSpeed * deltaTime);
			});
	if (key == GLFW_KEY_S && action == GLFW_RELEASE)
		scene.keys[BTNS].BindUpEvent([]() {
			});

	if (key == GLFW_KEY_A && action == GLFW_PRESS)
		scene.keys[BTNA].BindKeepEvent([]() {
		MainCamera::GetInstance()->LRMove(MainCamera::GetInstance()->cameraSpeed * deltaTime);
			});
	if (key == GLFW_KEY_A && action == GLFW_RELEASE)
		scene.keys[BTNA].BindUpEvent([]() {
			});

	if (key == GLFW_KEY_D && action == GLFW_PRESS)
		scene.keys[BTND].BindKeepEvent([]() {
		MainCamera::GetInstance()->LRMove(-MainCamera::GetInstance()->cameraSpeed * deltaTime);
			});
	if (key == GLFW_KEY_D && action == GLFW_RELEASE)
		scene.keys[BTND].BindUpEvent([]() {
			});

	if (key == GLFW_KEY_1 && action == GLFW_PRESS)
		scene.keys[BTN1].BindDownEvent([]() {
		scene.drawMode.isLine = !scene.drawMode.isLine;
			});
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
		scene.mouse.mouseRightDown = true;
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE)
		scene.mouse.mouseRightDown = false;
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	//cout << xpos << "  " << ypos << endl;
	if (scene.mouse.mouseRightDown)
	{
		float disx = scene.mouse.cursorPrePos.x - xpos;
		float disy = scene.mouse.cursorPrePos.y - ypos;
		MainCamera::GetInstance()->LRRotate(disx * deltaTime);
		MainCamera::GetInstance()->UDRotate(disy * deltaTime);
	}

	/*if (scene.mouse.isCatchPoint && scene.mouse.mouseLeftDown)
	{
		float disx = scene.mouse.cursorPrePos.x - xpos;
		float disy = scene.mouse.cursorPrePos.y - ypos;

		float temp = deltaTime * 2;
		translate(scene.ObjArray()[scene.mouse.catchObjIndex].World, vec3(disx * temp, disy * temp, 0.0));
	}*/

	scene.mouse.cursorNowPos = scene.mouse.cursorPrePos - vec2(xpos, ypos);

	scene.mouse.cursorPrePos.x = xpos;
	scene.mouse.cursorPrePos.y = ypos;
}