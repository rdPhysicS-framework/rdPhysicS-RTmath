#ifndef __RT_MATH_OPERATION_H__
#define __RT_MATH_OPERATION_H__

#include "Util.h"

namespace RT
{
	class Math
	{
	public:
		static inline double ToRadians(double degrees)
		{
			return degrees * RT_TO_RAD;
		}

		static inline float ToRadians(float degrees)
		{
			return static_cast<float>(degrees * RT_TO_RAD);
		}

		static inline double ToDegrees(double radians)
		{
			return radians * RT_TO_DEG;
		}

		static inline float ToDegrees(float radians)
		{
			return static_cast<float>(radians * RT_TO_DEG);
		}
	};
}
#endif //__RT_MATH_OPERATION_H__

