#include "RTvec4.h"
#include "RTmat3.h"
#include "RTmat4.h"
#include "RTquaternion.h"
#include <sstream>
#include <cmath>

using namespace RT;

Vec4f::Vec4f(float _x, float _y, float _z, float _w) :
		  x(_x), y(_y), z(_z), w(_w)
{
}

Vec4f::Vec4f(const Vec4f &other) :
		x(other.x), y(other.y),
		z(other.z), w(other.w)
{}

Vec4f::Vec4f(const Vec3f &vector3, float _w) :
		x(vector3.x), y(vector3.y),
		z(vector3.z), w(_w)
{}

Vec4f::Vec4f(const Vec2f &vector2, float _z, float _w) :
		x(vector2.x), y(vector2.y),
		z(_z), w(_w)
{}

Vec4f::Vec4f(float vector[SIZE_4]) :
		  x(vector[0]), y(vector[1]),
		  z(vector[2]), w(vector[3])
{
}

Vec4f::~Vec4f()
{}

Vec4f &Vec4f::Set(float _x, float _y, float _z, float _w)
{
	x = _x;
	y = _y;
	z = _z;
	w = _w;

	return (*this);
}

Vec4f &Vec4f::Set(const Vec4f &other)
{
	x = other.x;
	y = other.y;
	z = other.z;
	w = other.w;

	return (*this);
}

Vec4f &Vec4f::Set(const Vec3f &vector3f, float _w)
{
	x = vector3f.x;
	y = vector3f.y;
	z = vector3f.z;
	w = _w;

	return (*this);
}

Vec4f &Vec4f::Set(const Vec2f &other, float _z, float _w)
{
	x = other.x;
	y = other.y;
	z = _z;
	w = _w;

	return (*this);
}

Vec4f &Vec4f::Set(float vector[SIZE_4])
{
	x = vector[0];
	y = vector[1];
	z = vector[2];
	w = vector[3];

	return (*this);
}

Vec4f &Vec4f::operator=(const Vec4f &other)
{
	if (this != &other)
	{
		return Set(other);
	}

	return (*this);
}

Vec4f &Vec4f::operator=(const Vec3f &other)
{
	return Set(other, 1.0f);
}

Vec4f &Vec4f::operator=(const Vec2f &other)
{
	return Set(other, 0, 1);
}

Vec4f Vec4f::operator+(const Vec4f &other) const
{
	return Vec4f(x + other.x,
				  y + other.y,
				  z + other.z,
				  w + other.w);
}

Vec4f &Vec4f::operator+=(const Vec4f &other)
{
	x += other.x;
	y += other.y;
	z += other.z;
	w += other.w;

	return (*this);
}

Vec4f Vec4f::operator+(float vector[SIZE_4]) const
{
	return Vec4f(x + vector[0],
				  y + vector[1],
				  z + vector[2],
				  w + vector[3]);
}

Vec4f &Vec4f::operator+=(float vector[SIZE_4])
{
	x += vector[0];
	y += vector[1];
	z += vector[2];
	w += vector[3];

	return (*this);
}

Vec4f Vec4f::operator+(const Vec3f &other) const
{
	return Vec4f(x + other.x,
				  y + other.y,
				  z + other.z,
				  w); 
}

Vec4f &Vec4f::operator+=(const Vec3f &other)
{
	x += other.x;
	y += other.y;
	z += other.z;

	return (*this);
}

Vec4f Vec4f::operator+(const Vec2f &other) const
{
	return Vec4f(x + other.x,
				  y + other.y,
				  z, w);
}

Vec4f &Vec4f::operator+=(const Vec2f &other)
{
	x += other.x;
	y += other.y;

	return (*this);
}

Vec4f Vec4f::operator+(float scalar) const
{
	return Vec4f(x + scalar,
				  y + scalar,
				  z + scalar,
				  w + scalar);
}

Vec4f &Vec4f::operator+=(float scalar)
{
	x += scalar;
	y += scalar;
	z += scalar;
	w += scalar;

	return (*this);
}

Vec4f Vec4f::operator-(const Vec4f &other) const
{
	return Vec4f(x - other.x,
				  y - other.y,
				  z - other.z,
				  w - other.w);
}

Vec4f &Vec4f::operator-=(const Vec4f &other)
{
	x -= other.x;
	y -= other.y;
	z -= other.z;
	w -= other.w;

	return (*this);
}

Vec4f Vec4f::operator-(float vector[SIZE_4]) const
{
	return Vec4f(x - vector[0],
				  y - vector[1],
				  z - vector[2],
				  w - vector[3]);
}

Vec4f &Vec4f::operator-=(float vector[SIZE_4])
{
	x -= vector[0];
	y -= vector[1];
	z -= vector[2];
	w -= vector[3];

	return (*this);
}

Vec4f Vec4f::operator-(const Vec3f &other) const
{
	return Vec4f(x - other.x,
				  y - other.y,
				  z - other.z,
				  w);
}

Vec4f &Vec4f::operator-=(const Vec3f &other)
{
	x -= other.x;
	y -= other.y;
	z -= other.z;

	return (*this);
}

Vec4f Vec4f::operator-(const Vec2f &other) const
{
	return Vec4f(x - other.x,
				  y - other.y,
				  z, w);
}

Vec4f &Vec4f::operator-=(const Vec2f &other)
{
	x -= other.x;
	y -= other.y;

	return (*this);
}

Vec4f Vec4f::operator-(float scalar) const
{
	return Vec4f(x - scalar,
				  y - scalar,
				  z - scalar,
				  w - scalar);
}

Vec4f &Vec4f::operator-=(float scalar)
{
	x -= scalar;
	y -= scalar;
	z -= scalar;
	w -= scalar;

	return (*this);
}

Vec4f Vec4f::operator*(const Vec4f &other) const
{
	return Vec4f(x * other.x,
				  y * other.y,
				  z * other.z,
				  w * other.w);
}

Vec4f &Vec4f::operator*=(const Vec4f &other)
{
	x *= other.x;
	y *= other.y;
	z *= other.z;
	w *= other.w;

	return (*this);
}

Vec4f Vec4f::operator*(float vector[SIZE_4]) const
{
	return Vec4f(x * vector[0],
				  y * vector[1],
				  z * vector[2],
				  w * vector[3]);
}

Vec4f &Vec4f::operator*=(float vector[SIZE_4])
{
	x *= vector[0];
	y *= vector[1];
	z *= vector[2];
	w *= vector[3];

	return (*this);
}

Vec4f Vec4f::operator*(Mat4f &matrix4) const
{
	return (Vec4f(*this) *= matrix4);
}

Vec4f &Vec4f::operator*=(Mat4f &matrix4)
{
	float vet[SIZE_4] = { 0, 0, 0, 0 };

	for (int i = 0; i < SIZE_4; i++)
	{
		for (int j = 0; j < SIZE_4; j++)
		{
			vet[i] += (*this)[j] * matrix4[j][i];
		}
	}

	return Set(vet);
}

Vec4f Vec4f::operator*(float scalar) const
{
	return Vec4f(x * scalar,
				  y * scalar,
				  z * scalar,
				  w * scalar);;
}

Vec4f &Vec4f::operator*=(float scalar)
{
	x *= scalar;
	y *= scalar;
	z *= scalar;
	w *= scalar;

	return (*this);
}

Vec4f Vec4f::operator/(float scalar) const
{
	return Vec4f(x / scalar,
				  y / scalar,
				  z / scalar,
				  w / scalar);
}

Vec4f &Vec4f::operator/=(float scalar)
{
	x /= scalar;
	y /= scalar;
	z /= scalar;
	w /= scalar;

	return (*this);
}

Vec4f &Vec4f::Normalize()
{
	float value = 1 / Size();

	return (value == 0) ? (*this) : 
			Set(x*value, y*value, 
				z*value, w*value);
}

Vec4f &Vec4f::Lerp(const Vec4f & end, float t)
{
	return Set(x + (end.x - x) * t,
			   y + (end.y - y) * t,
			   z + (end.z - z) * t,
			   w + (end.w - w) * t);
}

Vec4f &Vec4f::Negate()
{
	return (*this *= -1);
}

Vec4f &RT::Vec4f::Mix(const Vec4f &other, float percentage)
{
	x = x * (1.0f - percentage) + other.x * percentage;
	y = y * (1.0f - percentage) + other.y * percentage;
	z = z * (1.0f - percentage) + other.z * percentage;
	w = w * (1.0f - percentage) + other.w * percentage;

	return (*this);
}

float Vec4f::GetAngle(const Vec4f &other) const
{
	float dot = Dot(other);
	float module = Size() * other.Size();

	return acosf(dot / module);
}

float Vec4f::Dot(const Vec4f &other) const
{
	return ((x * other.x) +
			(y * other.y) +
			(z * other.z) +
			(w * other.w));
}

float Vec4f::Size() const
{
	return sqrtf((x * x) +
				 (y * y) +
				 (z * z) +
				 (w * w));
}

float Vec4f::SizeSQR() const
{
	return	(x * x) +
			(y * y) +
			(z * z) +
			(w * w);
}

std::ostream &RT::operator<<(std::ostream &out, const Vec4f &vec4)
{
	out << "(" << vec4.x << ", " << vec4.y << ", " 
		       << vec4.z << ", " << vec4.w << ")" << std::endl;

	return out;
}

Vec4f RT::operator*(float value, const Vec4f &vec)
{
	return (vec * value);
}

Vec4f RT::operator+(float value, const Vec4f & vec)
{
	return (vec + value);
}

Vec4f RT::operator-(float value, const Vec4f & vec)
{
	return (vec - value);
}
