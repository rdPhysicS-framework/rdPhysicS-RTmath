#include "RTvec3.h"
#include "RTvec4.h"
#include "RTmat3.h"
#include "RTmathOp.h"
#include <iostream>
#include <cmath>
#include <sstream>

using namespace RT;

Vec3f::Vec3f(float _x, float _y, float _z)  : 
		  x(_x), y(_y), z(_z)
 {}

Vec3f::Vec3f(const Vec3f &other) :
	x(other.x), y(other.y), z(other.z)
{}

Vec3f::Vec3f(float vector[SIZE_3]) :
		  x(vector[0]), 
		  y(vector[1]), 
		  z(vector[2])
{}

Vec3f::Vec3f(const Vec2f &other, float _z) : 
	x(other.x), y(other.y), z(_z)
{}

Vec3f &Vec3f::Set(float _x, float _y, float _z)
{
	x = _x;
	y = _y;
	z = _z;

	return (*this);
}

Vec3f &Vec3f::Set(const Vec3f &other)
{
	x = other.x;
	y = other.y;
	z = other.z;

	return (*this);
}

Vec3f &Vec3f::Set(float vector[SIZE_3])
{
	x = vector[0];
	y = vector[1];
	z = vector[2];

	return (*this);
}

Vec3f &Vec3f::Set(const Vec2f &other, float _z)
{
	x = other.x;
	y = other.y;
	z = _z;

	return (*this);
}

Vec3f Vec3f::operator+(const Vec3f &other) const
{
	return Vec3f(x + other.x,
				  y + other.y, 
				  z + other.z);
}

Vec3f &Vec3f::operator+=(const Vec3f &other)
{
	x += other.x;
	y += other.y;
	z += other.z;
	
	return (*this);
}

Vec3f Vec3f::operator+(const Vec2f &other) const
{
	return Vec3f(x + other.x,
				  y + other.y, z);
}

Vec3f &Vec3f::operator+=(const Vec2f &other)
{
	x += other.x;
	y += other.y;

	return (*this);
}

Vec3f Vec3f::operator+(float value) const
{
	return Vec3f(x + value,
				  y + value,
				  z + value);
}

Vec3f &Vec3f::operator+=(float value) 
{
	x += value;
	y += value;
	z += value;

	return (*this);
}

Vec3f Vec3f::operator-(const Vec3f &other) const
{
	return Vec3f(x - other.x,
				  y - other.y,
				  z - other.z);
}

Vec3f &Vec3f::operator-=(const Vec3f &other)
{
	x -= other.x;
	y -= other.y;
	z -= other.z;
	
	return (*this);
}

Vec3f Vec3f::operator-(const Vec2f &other) const
{
	return Vec3f(x - other.x,
				  y - other.y,
				  z           );
}

Vec3f &Vec3f::operator-=(const Vec2f &other)
{
	this->x -= other.x;
	this->y -= other.y;

	return (*this);
}

Vec3f Vec3f::operator-(float value) const
{
	return Vec3f(x - value,
				  y - value,
				  z - value);
}

Vec3f &Vec3f::operator-=(float value)
{
	x -= value;
	y -= value;
	z -= value;

	return (*this);
}

Vec3f Vec3f::operator*(const Vec3f &other) const
{
	return Vec3f(x * other.x,
				  y * other.y, 
				  z * other.z);
}

Vec3f &Vec3f::operator*=(const Vec3f &other)
{
	x *= other.x;
	y *= other.y;
	z *= other.z;
	
	return (*this);
}

Vec3f Vec3f::operator*(Mat3f &matrix) const
{
	return (Vec3f(*this) *= matrix);
}

Vec3f &Vec3f::operator*=(Mat3f &matrix)
{
	float vec[SIZE_3] = { 0, 0, 0 };

	for (int i = 0, count = 0; i < SIZE_3; i++, count++)
	{
		for (int j = 0; j < SIZE_3; j++)
		{
			vec[count] += (*this)[j] * matrix[j][i];
		}
	}

	return Set(vec);
}

Vec3f Vec3f::operator*(float scalar) const
{
	return Vec3f(x * scalar,
				  y * scalar,
				  z * scalar);
}

Vec3f &Vec3f::operator*=(float scalar)
{
	x *= scalar;
	y *= scalar;
	z *= scalar;
	
	return (*this);
}

Vec3f Vec3f::operator/(float scalar) const
{
	return (scalar != 0)?
			Vec3f(x / scalar, 
				   y / scalar, 
				   z / scalar) :
		    throw std::exception("It is not possible to divide for zero\n");
}

Vec3f &Vec3f::operator/=(float scalar)
{
	if (scalar != 0)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;
	}
	else
		throw std::exception("It is not possible to divide for zero\n");

	return (*this);
}

Vec3f &Vec3f::operator=(const Vec3f &other)
{
	if (this != &other)
	{
		return Set(other);
	}

	return (*this);
}

Vec3f &Vec3f::operator=(const Vec2f &other)
{
	return Set(other);
}

float Vec3f::Dot(const Vec3f &other) const
{
	return (x * other.x) +
		   (y * other.y) + 
		   (z * other.z);

}

Vec3f &Vec3f::Cross(const Vec3f &other)
{
	return Set((y * other.z) - (z * other.y),
			   (z * other.x) - (x * other.z),
			   (x * other.y) - (y * other.x));
}

Vec3f &Vec3f::Normalize()
{
	float value = 1.0f / Size();

	return (value == 0)? (*this) : 
		   Set(x*value, y*value, z*value);
}

Vec3f &Vec3f::Mix(const Vec3f &other, float percentage)
{
	x = x * (1.0f - percentage) + other.x * percentage;
	y = y * (1.0f - percentage) + other.y * percentage;
	z = z * (1.0f - percentage) + other.z * percentage;

	return (*this);
}

float Vec3f::Size() const
{
	return sqrtf(x*x + y*y + z*z);
}

float Vec3f::SizeSQR() const
{
	return (x * x) + 
		   (y * y) + 
		   (z * z);
}

Vec3f &Vec3f::RotateX(float radians)
{
	float Sin = sinf(radians);
	float Cos = cosf(radians);

	float y = (y * Cos) - (z * Sin);
	float z = (z * Cos) + (y * Sin);

	y = y;
	z = z;

	return (*this);
}

Vec3f &Vec3f::RotateY(float radians)
{
	float Sin = sinf(radians);
	float Cos = cosf(radians);

	float x = (x * Cos) + (z * Sin);
	float z = (z * Cos) - (x * Sin);

	x = x;
	z = z;

	return (*this);
}

Vec3f &Vec3f::RotateZ(float radians)
{
	float Sin = sinf(radians);
	float Cos = cosf(radians);

	float x = (x * Cos) - (y * Sin);
	float y = (y * Cos) + (x * Sin);

	x = x;
	y = y;

	return (*this);
}

Vec3f &RT::Vec3f::Rotate(const Vec3f &axis, float radians)
{
	float Sin = sinf(radians);
	float Cos = cosf(radians);
	float k = 1.0f - Cos;

	x = x * (Cos + k * axis.x * axis.x) +
		y * (k * axis.x * axis.y - Sin * axis.z) +
		z * (k * axis.x * axis.z + Sin * axis.x);

	y = x * (k * axis.x * axis.y + Sin * axis.z) +
		y * (Cos + k * axis.y * axis.y) +
		z * (k * axis.y * axis.z - Sin * axis.x);

	z = x * (k * axis.x * axis.z - Sin * axis.y) +
		y * (k * axis.y * axis.z + Sin * axis.x) +
		z * (Cos + k * axis.z * axis.z);

	return(*this);
}

Vec3f &Vec3f::Reflect(const Vec3f &normal)
{
	float dot = Dot(normal);
	x = x - 2.0f * dot * normal.x;
	y = y - 2.0f * dot * normal.y;
	z = z - 2.0f * dot * normal.z;

	return (*this);
}

Vec3f &Vec3f::Refract(const Vec3f &normal, float index)
{
	float dot = Dot(normal);
	float k = 1.0f - pow(index, 2) * (1.0f - pow(dot, 2));
	
	if (k >= 0.0f) 
	{
		float r = sqrtf(k);
		x = index * x - normal.x * (index * dot + r);
		y = index * y - normal.y * (index * dot + r);
		z = index * z - normal.z * (index * dot + r);
	}

	return (*this);
}

Vec3f &Vec3f::CartesianCoordinates(float ray, float angleZenite, float angleAzimute)
{
	x = ray * cosf(angleZenite) * cosf(angleAzimute);
	y = ray * sinf(angleZenite);
	z = ray * sinf(angleAzimute) * cosf(angleZenite);

	return (*this);
}

Vec3f &Vec3f::PolarCoordinate()
{
	float ray = Size();
	float zenite = Math::ToDegrees(atan2f(z, x));
	float azimute = Math::ToDegrees(atan2f(Vec2f(x, y).Size(), y));

	return Set(ray, zenite, azimute);
}

Vec3f &Vec3f::Lerp(const Vec3f &other, float t)
{
	x = x + (other.x - x) * t;
	y = y + (other.y - y) * t;
	z = z + (other.z - z) * t;

	return (*this);
}

std::ostream &RT::operator<<(std::ostream &out, const Vec3f &vec3)
{
	out << "(" << vec3.x << ", " << vec3.y << ", "
			   << vec3.z <<  ")" << std::endl;

	return out;
}

Vec3f RT::operator*(float value, const Vec3f &vec)
{
	return (vec * value);
}

Vec3f RT::operator+(float value, const Vec3f &vec)
{
	return (vec + value);
}

Vec3f RT::operator-(float value, const Vec3f &vec)
{
	return (vec - value);
}
