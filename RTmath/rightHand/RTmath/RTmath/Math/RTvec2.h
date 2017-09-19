#ifndef __RT_VEC2_H__
#define __RT_VEC2_H__

#include "Util.h"
#include <ostream>
#include <sstream>

namespace RT
{
	class Vec2f
	{
	public:
		float x;
		float y;

	public:
		explicit Vec2f(float _x = 0, float _y = 0);
		Vec2f(const Vec2f &other);
		explicit Vec2f(float other[SIZE_2]);

		Vec2f &Set(float _x, float _y);
		Vec2f &Set(const Vec2f &other);
		Vec2f &Set(float other[SIZE_2]);

		inline float Vec2f::operator[](const int i) const
		{
			if (i == 0) return x;
			if (i == 1) return y;

			std::stringstream s;
			s << "Index out of bounds " << i << "size 2.\n";
			throw std::out_of_range(s.str());
		}

		inline float Vec2f::operator+(const int i) const
		{
			if (i == 0) return x;
			if (i == 1) return y;

			std::stringstream s;
			s << "Index out of bounds " << i << "size 2.\n";
			throw std::out_of_range(s.str());
		}

		inline float &Vec2f::operator[](const int i)
		{
			if (i == 0) return x;
			if (i == 1) return y;

			std::stringstream s;
			s << "Index out of bounds " << i << "size 2.\n";
			throw std::out_of_range(s.str());
		}

		inline float &Vec2f::operator+(const int i)
		{
			if (i == 0) return x;
			if (i == 1) return y;

			std::stringstream s;
			s << "Index out of bounds " << i << "size 2.\n";
			throw std::out_of_range(s.str());
		}

		Vec2f operator+(const Vec2f &other) const;
		Vec2f &operator+=(const Vec2f &other);
		Vec2f operator+(float value) const;
		Vec2f &operator+=(float value);
		Vec2f operator-(const Vec2f &other) const;
		Vec2f &operator-=(const Vec2f &other);
		Vec2f operator-(float value) const;
		Vec2f &operator-=(float value);
		Vec2f operator*(const Vec2f &other) const;
		Vec2f &operator*=(const Vec2f &other);
		Vec2f operator*(float value) const;
		Vec2f &operator*=(float value);
		Vec2f operator/(float value) const;
		Vec2f &operator/=(float value);
		Vec2f &operator=(const Vec2f &other);

		inline bool Vec2f::operator==(const Vec2f &other) const
		{
			return x == other.x &&
				   y == other.y;
		}

		inline bool Vec2f::operator!=(const Vec2f &other) const
		{
			return x != other.x ||
				   y != other.y;
		}

		Vec2f &Normalize();
		Vec2f &Rotate(float radius);
		Vec2f &Refract(const Vec2f &normal, float index);
		Vec2f &Reflect(const Vec2f &normal);
		Vec2f &Lerp(const Vec2f &end, float t);

		float GetAngle() const;
		float Dot(const Vec2f &other) const;
		float Size() const;
		float SizeSQR() const;
	};

	std::ostream &operator<<(std::ostream &out, const Vec2f &vec2);
	Vec2f operator*(float value, const Vec2f &vec);
	Vec2f operator+(float value, const Vec2f &vec);
	Vec2f operator-(float value, const Vec2f &vec);

	namespace vc2
	{
		inline float Dot(const Vec2f &vec1, const Vec2f &vec2)
		{
			return vec1.Dot(vec2);
		}

		inline Vec2f Lerp(const Vec2f &start, const Vec2f &end, float t)
		{
			return Vec2f(start).Lerp(end, t);
		}

		inline Vec2f Normalize(const Vec2f &vec)
		{
			return Vec2f(vec).Normalize();
		}

		inline Vec2f Refract(const Vec2f &vec1, const Vec2f &vec2, float index)
		{
			return Vec2f(vec1).Refract(vec2, index);
		}

		inline Vec2f Reflect(const Vec2f &vec, const Vec2f &normal)
		{
			return Vec2f(vec).Reflect(normal);
		}

		inline Vec2f Rotate(const Vec2f &vec, float angle)
		{
			return Vec2f(vec).Rotate(angle);
		}

		inline Vec2f Negate(const Vec2f &vec)
		{
			return (vec * -1);
		}

		inline Vec2f Perp(const Vec2f & vec)
		{
			return Vec2f(-vec.y, vec.x);
		}
	}
}
#endif //__RT_VEC2_H__