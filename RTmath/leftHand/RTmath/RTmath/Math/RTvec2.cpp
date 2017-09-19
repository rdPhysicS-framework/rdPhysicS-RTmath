#include "RTvec2.h"
#include <iostream>
#include <cmath>

using namespace RT;

Vec2f::Vec2f(float _x, float _y) : 
	x(_x), y(_y)
{}

Vec2f::Vec2f(const Vec2f &other) :
	x(other.x), y(other.y)
{}

Vec2f::Vec2f(float other[SIZE_2]) : 
	x(other[0]), y(other[1])
{}

Vec2f &Vec2f::Set(float _x, float _y)
{
	x = _x;
	y = _y;

	return (*this);
}

Vec2f &Vec2f::Set(const Vec2f &other)
{
	x = other.x;
	y = other.y;

	return (*this);
}

Vec2f &Vec2f::Set(float other[SIZE_2])
{
	x = other[0];
	y = other[1];

	return (*this);
}

Vec2f Vec2f::operator+(const Vec2f &other) const
{
	return Vec2f(x + other.x,
				  y + other.y);
}

Vec2f &Vec2f::operator+=(const Vec2f &other)
{
	x += other.x;
	y += other.y;

	return (*this);
}

Vec2f Vec2f::operator+(float value) const
{
	return Vec2f(x + value,
				  y + value);
}

Vec2f &Vec2f::operator+=(float value)
{
	x += value;
	y += value;

	return (*this);
}

Vec2f Vec2f::operator-(const Vec2f &other) const
{
	return Vec2f(x - other.x,
				  y - other.y);
}

Vec2f &Vec2f::operator-=(const Vec2f &other)
{
	x -= other.x;
	y -= other.y;

	return (*this);
}

Vec2f Vec2f::operator-(float value) const
{
	return Vec2f(x - value,
				  y - value);
}

Vec2f &Vec2f::operator-=(float value)
{
	x -= value;
	y -= value;

	return (*this);
}

Vec2f Vec2f::operator*(const Vec2f &other) const
{
	return Vec2f(x * other.x,
				  y * other.y);
}

Vec2f &Vec2f::operator*=(const Vec2f &other)
{
	x *= other.x;
	y *= other.y;

	return (*this);
}

Vec2f Vec2f::operator*(float value) const
{
	return Vec2f(x * value,
				  y * value);
}

Vec2f &Vec2f::operator*=(float value)
{
	x *= value;
	y *= value;

	return (*this);
}

Vec2f Vec2f::operator/(float value) const
{
	return Vec2f(x / value, y / value);
}

Vec2f &Vec2f::operator/=(float value)
{
	x /= value;
	y /= value;

	return (*this);
}

Vec2f &Vec2f::operator=(const Vec2f &other)
{
	if (*this != other)
	{
		x = other.x;
		y = other.y;
	}

	return (*this);
}

Vec2f &Vec2f::Normalize()
{
	float value = 1 / Size();
	
	return (value == 1E-7f)? (*this) : Set(x*value, y*value);
}

Vec2f &Vec2f::Rotate(float radius)
{
	float Sin = sinf(radius);
	float Cos = cosf(radius);

	float x = x * Cos - y * Sin;
	float y = x * Sin + y * Cos;

	return Set(x, y);
}

Vec2f &Vec2f::Refract(const Vec2f &normal, float index)
{
	float dot = Dot(normal);
	float k = 1.0f - powf(index, 2) * (1.0f * powf(dot, 2));

	if (k >= 0.0f)
	{
		float r = sqrt(k);
		x = index * x - normal.x * (index * dot * r);
		y = index * y - normal.y * (index * dot * r);
	}

	return (*this);
}

Vec2f &Vec2f::Reflect(const Vec2f &normal)
{
	float dot = Dot(normal);
	return Set(x - normal.x * 2.0f * dot,
			   y - normal.y * 2.0f * dot);
}

Vec2f &Vec2f::Lerp(const Vec2f &end, float t)
{
	return Set(x + t * (end.x - x),
			   y + t * (end.y - y));
}

float Vec2f::GetAngle() const
{
	return atan2f(y, x);
}

float Vec2f::Dot(const Vec2f &other) const
{
	return (x * other.x) +
		   (y * other.y);
}

float Vec2f::Size() const
{
	return sqrtf((x * x) + (y * y));
}

float Vec2f::SizeSQR() const
{
	return ((x * x) + (y * y));
}

std::ostream &RT::operator<<(std::ostream &out, const Vec2f &vec2)
{
	out << "(" << vec2.x << ", " << vec2.y << ")" << std::endl;

	return out;
}

Vec2f RT::operator*(float value, const Vec2f &vec)
{
	return (vec * value);
}

Vec2f RT::operator+(float value, const Vec2f &vec)
{
	return (vec + value);
}

Vec2f RT::operator-(float value, const Vec2f &vec)
{
	return (vec - value);
}

