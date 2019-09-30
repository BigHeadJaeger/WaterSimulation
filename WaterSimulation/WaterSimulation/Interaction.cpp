#include"Interaction.h"
void Key::EventDown()
{
	//执行完了之后解除绑定
	outerFunction();
	UnBind();
}

void Key::EventKeep()
{
	//isDown = true;
	outerFunction();					//不断执行
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
