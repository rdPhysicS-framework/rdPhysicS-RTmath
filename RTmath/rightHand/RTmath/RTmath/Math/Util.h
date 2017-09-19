#ifndef __UTIL_H__
#define __UTIL_H__

#include<cmath>

#define SIZE_2 2
#define SIZE_3 3
#define SIZE_4 4

#ifndef RT_ANGLE
	#define RT_ANGLE 180
#endif //RT_ANGLE

#ifndef RT_PI
	#define RT_PI 3.14159265358979323846
#endif //RT_PI

#ifndef RT_2PI
	#define RT_2PI (RT_PI * 2)
#endif //RT_PI

#ifndef RT_TO_DEG
	#define RT_TO_DEG (RT_ANGLE/RT_PI)
#endif //RT_TO_DEG

#ifndef RT_TO_RAD
	#define RT_TO_RAD  (RT_PI/RT_ANGLE)
#endif //RT_TO_RAD

#endif //__UTIL_H__
