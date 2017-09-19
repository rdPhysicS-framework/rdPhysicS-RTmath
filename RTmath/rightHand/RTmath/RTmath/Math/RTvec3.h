#ifndef __RT_VEC3_H__
#define __RT_VEC3_H__

#include "RTvec2.h"

namespace RT
{
	class Mat3f;

	class Vec3f
	{
	public:
		float x, y, z;

	public:
		explicit Vec3f(float _x = 0, float _y = 0, float _z = 0);
		Vec3f(const Vec3f &other);
		explicit Vec3f(float vector[SIZE_3]);
		explicit Vec3f(const Vec2f &other, float _z = 0);

		inline Vec2f XY() const { return Vec2f(x, y); }
		inline Vec2f XZ() const { return Vec2f(x, z); }
		inline Vec2f YX() const { return Vec2f(y, x); }
		inline Vec2f YZ() const { return Vec2f(y, z); }
		inline Vec2f ZX() const { return Vec2f(z, x); }
		inline Vec2f ZY() const { return Vec2f(z, y); }

		Vec3f &Set(float _x, float _y, float _z);
		Vec3f &Set(const Vec3f &other);
		Vec3f &Set(float vector[SIZE_3]);
		Vec3f &Set(const Vec2f &other, float _z = 0);

		Vec3f operator+(const Vec3f &other) const;
		Vec3f &operator+=(const Vec3f &other);
		Vec3f operator+(const Vec2f &other) const;
		Vec3f &operator+=(const Vec2f &other);
		Vec3f operator+(float value) const;
		Vec3f &operator+=(float value);
		Vec3f operator-(const Vec3f &other) const;
		Vec3f &operator-=(const Vec3f &other);
		Vec3f operator-(const Vec2f &other) const;
		Vec3f &operator-=(const Vec2f &other);
		Vec3f operator-(float value) const;
		Vec3f &operator-=(float value);
		Vec3f operator*(const Vec3f &other) const;
		Vec3f &operator*=(const Vec3f &other);
		Vec3f operator*(Mat3f &matrix) const;
		Vec3f &operator*=(Mat3f &matrix);
		Vec3f operator*(float scalar) const;
		Vec3f &operator*=(float scalar);
		Vec3f operator/(float scalar) const;
		Vec3f &operator/=(float scalar);

		float operator[](const int index) const
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;

			std::stringstream s;
			s << "Index out of bounds " << index << "size 3.\n";
			throw std::out_of_range(s.str());
		}

		float operator+(const int index) const
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;

			std::stringstream s;
			s << "Index out of bounds " << index << "size 3.\n";
			throw std::out_of_range(s.str());
		}

		float &operator[](const int index)
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;

			std::stringstream s;
			s << "Index out of bounds " << index << "size 3.\n";
			throw std::out_of_range(s.str());
		}

		float &operator+(const int index)
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;

			std::stringstream s;
			s << "Index out of bounds " << index << "size 3.\n";
			throw std::out_of_range(s.str());
		}

		Vec3f &operator=(const Vec3f &other);
		Vec3f &operator=(const Vec2f &other);

		inline bool Vec3f::operator==(const Vec3f &other) const
		{
			return (x == other.x &&
				    y == other.y &&
				    z == other.z);
		}

		inline bool Vec3f::operator!=(const Vec3f &other) const
		{
			return (x != other.x ||
				    y != other.y ||
				    z != other.z);
		}

		float Dot(const Vec3f &other) const;
		Vec3f &Cross(const Vec3f &other);
		Vec3f &Normalize();
		Vec3f &Mix(const Vec3f &other, float percentage);

		float Size() const;

		float SizeSQR() const;

		Vec3f &RotateX(float radians);
		Vec3f &RotateY(float radians);
		Vec3f &RotateZ(float radians);
		Vec3f &Rotate(const Vec3f &axis, float radians);

		Vec3f &Reflect(const Vec3f &normal);
		Vec3f &Refract(const Vec3f &normal, float index);
		Vec3f &CartesianCoordinates(float ray, float angleZenite, float angleAzimute);
		Vec3f &PolarCoordinate();

		Vec3f &Lerp(const Vec3f &other, float t);
	};

	std::ostream &operator<<(std::ostream &out, const Vec3f &vec3);
	Vec3f operator*(float value, const Vec3f &vec);
	Vec3f operator+(float value, const Vec3f &vec);
	Vec3f operator-(float value, const Vec3f &vec);

	namespace vc3
	{
		inline Vec3f CartesianCoordinates(float ray, float angleZenite, float angleAzimute)
		{
			return Vec3f().CartesianCoordinates(ray, angleZenite, angleAzimute);
		}

		inline Vec3f Cross(const Vec3f &vec1, const Vec3f &vec2)
		{
			return Vec3f((vec1.y * vec2.z) - (vec1.z * vec2.y),
						 (vec1.z * vec2.x) - (vec1.x * vec2.z),
						 (vec1.x * vec2.y) - (vec1.y * vec2.x));
		}

		inline float Dot(const Vec3f &vec1, const Vec3f &vec2)
		{
			return (vec1.x * vec2.x) +
				   (vec1.y * vec2.y) +
				   (vec1.z * vec2.z);
		}

		inline Vec3f Lerp(const Vec3f &start, const Vec3f &end, float t)
		{
			return Vec3f(start).Lerp(end, t);
		}

		inline Vec3f Mix(const Vec3f &vec1, const Vec3f &vec2, float percentage)
		{
			return Vec3f(vec1).Mix(vec2, percentage);
		}

		inline Vec3f Negate(Vec3f &vec)
		{
			return (vec * -1);
		}

		inline Vec3f Normalize(const Vec3f &vec)
		{
			return Vec3f(vec).Normalize();
		}

		inline Vec3f PolarCoordinate(const Vec3f &vec)
		{
			return Vec3f(vec).PolarCoordinate();
		}

		inline Vec3f Reflect(const Vec3f &vec, const Vec3f &normal)
		{
			return Vec3f(vec).Reflect(normal);
		}

		inline Vec3f Refract(const Vec3f &vec, const Vec3f &normal, float index)
		{
			return Vec3f(vec).Refract(normal, index);
		}

		inline Vec3f RotateX(const Vec3f &vec, float radius)
		{
			return Vec3f(vec).RotateX(radius);
		}

		inline Vec3f RotateY(const Vec3f &vec, float radius)
		{
			return Vec3f(vec).RotateY(radius);
		}

		inline Vec3f Rotate(const Vec3f &vec, const Vec3f &axis, float radians)
		{
			return Vec3f(vec).Rotate(axis, radians);
		}

		inline Vec3f RotateZ(const Vec3f &vec, float radius)
		{
			return Vec3f(vec).RotateZ(radius);
		}
	}
}
#endif //__RT_VEC3_H__
