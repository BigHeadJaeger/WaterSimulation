#include"Interaction.h"
void Key::EventDown()
{
	//ִ������֮������
	outerFunction();
	UnBind();
}

void Key::EventKeep()
{
	//isDown = true;
	outerFunction();					//����ִ��
}

void Key::EventUP()
{
	outerFunction();
	UnBind();
}

void Key::Execute()
{
	if (eventFunctionPtr != NULL)
	{
		eventFunctionPtr();
	}
}
