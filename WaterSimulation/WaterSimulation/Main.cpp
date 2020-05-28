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
	//创建窗口
	GLFWwindow* window;

	//初始化glfw
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(WIDTH, HEIGHT, "WaterSimulation", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* 该函数指定在当前调用线程上的窗口渲染环境为OpenGL或OpenGL ES环境 */
	glfwMakeContextCurrent(window);

	//指定事件
	glfwSetKeyCallback(window, key_callback);			//键盘事件
	//glfwSetInputMode(window, GLFW_STICKY_KEYS, 1);

	glfwSetMouseButtonCallback(window, mouse_button_callback);		//鼠标按键事件

	glfwSetCursorPosCallback(window, cursor_position_callback);		//鼠标指针事件


	//创建场景

	scene.Init();


	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		//计算一帧的间隔
		float currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;



		//场景更新和绘制
		scene.Update(deltaTime);
		scene.Draw();

		//交换双缓冲
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
	//绑定keep事件后一定要绑定up事件
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
		MainCamera::GetInstance()->LRRotate(disx * deltaTime * 0.5);
		MainCamera::GetInstance()->UDRotate(disy * deltaTime * 0.5);
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