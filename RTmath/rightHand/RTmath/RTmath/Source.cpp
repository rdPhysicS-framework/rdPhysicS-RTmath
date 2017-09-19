#include <iostream>
#include "RTmath.h"


int main()
{
	RT::Mat3f m(1, 2, 3, 
		10, 0, 4, 
		12, 0, 100);
	RT::Mat4f m2;
	m2.Identity();
	m2.Translate(1, 10, 2);

	std::cout << m2 << std::endl;

	float r = m[0][8];
	//std::cout << r << std::endl;
	std::cout << RT::mt4::Perspective(RT::Math::ToRadians(60.0f),
									  1920.0f/1080.0f, 0.1f, 1000.0f)<< std::endl;


	system("pause");
	return 0;
}