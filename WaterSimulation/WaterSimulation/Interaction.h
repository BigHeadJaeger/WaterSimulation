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

//����״̬����Ϣ
struct Mouse
{
	vec2 cursorPrePos;			//��¼����ǰһλ��
	vec2 cursorNowPos;			//��¼���ڵ�λ��
	bool mouseLeftDown;
	bool mouseRightDown;

	bool isSelect;				//��ʾ��ǰ���Ѿ�ѡ��Ķ���
	bool isCatchPoint;			//��ʾ��ǰ�����϶�����

	int catchObjIndex;			//��ǰ�϶��������index



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


	//�˺�����update�в���ִ�д�������function
	void Execute();

	//ָ��keyִ��down�¼�������ָ��һ����Ҫ�ⲿ����
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


