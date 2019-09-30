#pragma once
#include<functional>
#include<glm.hpp>
using namespace glm;
using namespace std;

enum KEYNAME
{
	BTNW,
	BTNA,
	BTNS,
	BTND,
	BTN1,
	NONE,
};

//鼠标的状态和信息
struct Mouse
{
	vec2 cursorPrePos;			//记录光标的前一位置
	vec2 cursorNowPos;			//记录现在的位置
	bool mouseLeftDown;
	bool mouseRightDown;

	bool isSelect;				//表示当前有已经选择的东西
	bool isCatchPoint;			//表示当前正在拖动顶点

	int catchObjIndex;			//当前拖动的物体的index



	Mouse() :cursorPrePos(vec2(0, 0)), mouseLeftDown(false), mouseRightDown(false)
	{
		isCatchPoint = false;
		isSelect = false;
	}
};

class Key
{
public:
	KEYNAME keyName;
	//bool isDown;
	function<void()> eventFunctionPtr;
	function<void()> outerFunction;
private:
	void EventKeep();
	void EventUP();
	void EventDown();

public:
	Key()
	{
		keyName = NONE;
		eventFunctionPtr = NULL;
		outerFunction = NULL;
	}
	Key(KEYNAME name)
	{
		keyName = name;
		eventFunctionPtr = NULL;
		outerFunction = NULL;
	}


	//此函数在update中不断执行传进来的function
	void Execute();

	//指定key执行down事件，并且指定一个需要外部函数
	template<typename F>
	void BindDownEvent(F f)
	{
		eventFunctionPtr = bind(&Key::EventDown, this);
		outerFunction = bind(f);
	}

	template<typename F>
	void BindKeepEvent(F f)
	{
		eventFunctionPtr = bind(&Key::EventKeep, this);
		outerFunction = bind(f);
	}

	template<typename F>
	void BindUpEvent(F f)
	{
		eventFunctionPtr = bind(&Key::EventUP, this);
		outerFunction = bind(f);
	}

	void UnBind()
	{
		eventFunctionPtr = NULL;
		outerFunction = NULL;
	}

	
};


