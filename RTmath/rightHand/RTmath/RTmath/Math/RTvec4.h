#ifndef __RT_VEC4_H__
#define __RT_VEC4_H__

#include "Util.h"
#include <ostream>
#include <sstream>

#include "RTvec3.h"
#include "RTvec2.h"

namespace RT
{
	class Mat3f;
	class Mat4f;
	class Quaternionf;

	class Vec4f
	{
	public:
		float x;
		float y;
		float z;
		float w;

	public:
		explicit Vec4f(float _x = 0, float _y = 0, float _z = 0, float _w = 1);
		Vec4f(const Vec4f &other);
		explicit Vec4f(const Vec3f &vector3, float w = 1);
		explicit Vec4f(const Vec2f &vector2, float z = 0, float w = 1);
		explicit Vec4f(float vector[SIZE_4]);
		~Vec4f();

		inline Vec3f XYZ() const { return Vec3f(x, y, z); }
		inline Vec3f YXZ() const { return Vec3f(y, x, z); }
		inline Vec3f YZX() const { return Vec3f(y, z, x); }
		inline Vec3f ZXY() const { return Vec3f(z, x, y); }
		inline Vec3f ZYX() const { return Vec3f(z, y, x); }
		inline Vec2f XY() const { return Vec2f(x, y); }
		inline Vec2f XZ() const { return Vec2f(x, z); }
		inline Vec2f Xw() const { return Vec2f(x, w); }
		inline Vec2f YX() const { return Vec2f(y, x); }
		inline Vec2f YZ() const { return Vec2f(y, z); }
		inline Vec2f YW() const { return Vec2f(y, w); }
		inline Vec2f ZX() const { return Vec2f(z, x); }
		inline Vec2f ZY() const { return Vec2f(z, y); }
		inline Vec2f ZW() const { return Vec2f(z, w); }
		inline Vec2f WX() const { return Vec2f(w, x); }
		inline Vec2f WY() const { return Vec2f(w, y); }
		inline Vec2f WZ() const { return Vec2f(w, z); }

		Vec4f &Set(float _x, float _y, float _z, float _w);
		Vec4f &Set(const Vec4f &other);
		Vec4f &Set(const Vec3f &vector3, float _w);
		Vec4f &Set(const Vec2f &vector2, float _z, float _w);
		Vec4f &Set(float vector[SIZE_4]);

		Vec4f &operator=(const Vec4f &other);
		Vec4f &operator=(const Vec3f &other);
		Vec4f &operator=(const Vec2f &other);

		inline bool Vec4f::operator==(const Vec4f &other) const
		{
			return (x == other.x &&
				    y == other.y &&
				    z == other.z &&
				    w == other.w);
		}

		inline bool Vec4f::operator!=(const Vec4f &other) const
		{
			return (x != other.x ||
				    y != other.y ||
				    z != other.z ||
				    w != other.w);
		}

		inline float operator[](const int index) const
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;
			if (index == 3) return w;

			std::stringstream s;
			s << "Index out of bounds " << index << " size 4.\n";
			throw std::out_of_range(s.str());
		}

		inline float operator+(const int index) const
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;
			if (index == 3) return w;

			std::stringstream s;
			s << "Index out of bounds " << index << " size 4.\n";
			throw std::out_of_range(s.str());
		}

		inline float &operator[](const int index)
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;
			if (index == 3) return w;

			std::stringstream s;
			s << "Index out of bounds " << index << " size 4.\n";
			throw std::out_of_range(s.str());
		}

		inline float &operator+(const int index)
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;
			if (index == 3) return w;

			std::stringstream s;
			s << "Index out of bounds " << index << " size 4.\n";
			throw std::out_of_range(s.str());
		}

		Vec4f operator+(const Vec4f &other) const;
		Vec4f &operator+=(const Vec4f &other);
		Vec4f operator+(float vector[SIZE_4]) const;
		Vec4f &operator+=(float vector[SIZE_4]);
		Vec4f operator+(const Vec3f &other) const;
		Vec4f &operator+=(const Vec3f &other);
		Vec4f operator+(const Vec2f &other) const;
		Vec4f &operator+=(const Vec2f &other);
		Vec4f operator+(float scalar) const;
		Vec4f &operator+=(float scalar);

		Vec4f operator-(const Vec4f &other) const;
		Vec4f &operator-=(const Vec4f &other);
		Vec4f operator-(float vector[SIZE_4]) const;
		Vec4f &operator-=(float vector[SIZE_4]);
		Vec4f operator-(const Vec3f &other) const;
		Vec4f &operator-=(const Vec3f &other);
		Vec4f operator-(const Vec2f &other) const;
		Vec4f &operator-=(const Vec2f &other);
		Vec4f operator-(float scalar) const;
		Vec4f &operator-=(float scalar);

		Vec4f operator*(const Vec4f &other) const;
		Vec4f &operator*=(const Vec4f &other);
		Vec4f operator*(float vector[SIZE_4]) const;
		Vec4f &operator*=(float vector[SIZE_4]);
		Vec4f operator*(Mat4f &matrix) const;
		Vec4f &operator*=(Mat4f &matrix);
		Vec4f operator*(float scalar) const;
		Vec4f &operator*=(float scalar);

		Vec4f operator/(float scalar) const;
		Vec4f &operator/=(float scalar);

		Vec4f &Normalize();
		Vec4f &Rotate(const Quaternionf &quaternion);
		Vec4f &Lerp(const Vec4f &end, float t);
		Vec4f &Negate();
		Vec4f &Mix(const Vec4f &other, float percentage);

		float GetAngle(const Vec4f &other) const;
		float Dot(const Vec4f &other) const;
		float Size() const;
		float SizeSQR() const;
	};

	std::ostream &operator<<(std::ostream &out, const Vec4f &vec4);
	Vec4f operator*(float value, const Vec4f &vec);
	Vec4f operator+(float value, const Vec4f &vec);
	Vec4f operator-(float value, const Vec4f &vec);

	namespace vc4
	{
		inline float Dot(const Vec4f &vec1, const Vec4f &vec2)
		{
			return ((vec1.x * vec2.x) +
				(vec1.y * vec2.y) +
				(vec1.z * vec2.z) +
				(vec1.w * vec2.w));
		}

		inline Vec4f Lerp(const Vec4f &start, const Vec4f &end, float t)
		{
			return Vec4f(start).Lerp(end, t);
		}

		inline Vec4f Mix(const Vec4f &vec1, const Vec4f &vec2, float percentage)
		{
			return Vec4f(vec1).Mix(vec2, percentage);
		}

		inline Vec4f Negate(const Vec4f &vec)
		{
			return (vec * -1);
		}

		inline Vec4f Normalize(const Vec4f &vec)
		{
			return Vec4f(vec).Normalize();
		}
	}
}
#endif //__RT_VEC4_H__

