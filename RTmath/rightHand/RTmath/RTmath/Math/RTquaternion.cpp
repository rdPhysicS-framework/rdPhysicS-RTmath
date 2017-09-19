#include "RTquaternion.h"
//#include "RTmathOp.h"
#include "RTmat4.h"
#include <sstream>

using namespace RT;

Quaternionf::Quaternionf(float _x, float _y, float _z, float _w) :
			x(_x), y(_y), z(_z), w(_w)
{}

Quaternionf::Quaternionf(float vector[SIZE_4]) :
			 x(vector[0]), y(vector[1]),
			 z(vector[2]), w(vector[3])
{
}

Quaternionf::Quaternionf(const Vec4f &vec4) :
	x(vec4.x), y(vec4.y), z(vec4.z), w(vec4.w)
{}

Quaternionf::Quaternionf(float angle, const Vec3f &axis)
{
	RotationAxis(angle, axis);
}

Quaternionf::Quaternionf(float angle1, const Vec3f &axis1, 
						   float angle2, const Vec3f &axis2, 
						   float angle3, const Vec3f &axis3)
{
	RotationAxis(angle1, axis1,
				 angle2, axis2,
				 angle3, axis3);
}

Quaternionf::Quaternionf(const Quaternionf &other) :
	x(other.x), y(other.y), z(other.z), w(other.w)
{}

Quaternionf::~Quaternionf()
{}

Quaternionf &Quaternionf::RotationAxis(float angle, const Vec3f &axis)
{
	return RotationAxis(angle  , axis.x, 
					    axis.y , axis.z);
}

Quaternionf &Quaternionf::RotationAxis(float angle, float _x,
									     float _y,    float _z)
{
	float Sin = sinf(angle * 0.5f);

	float inv = 1.0f / sqrtf(_x * _x + 
							 _y * _y + 
							 _z * _z);

	x = _x * Sin * inv;
	y = _y * Sin * inv;
	z = _z * Sin * inv;
	w = cosf(angle * 0.5f);

	return (*this);
}

Quaternionf &Quaternionf::RotationAxis(float angle1, const Vec3f &axis1,
									     float angle2, const Vec3f &axis2, 
									     float angle3, const Vec3f &axis3)
{
	Quaternionf q1(angle1, axis1);
	Quaternionf q2(angle2, axis2);
	Quaternionf q3(angle3, axis3);

	*this = (q1 * q2 * q3);

	return (*this);
}

inline Quaternionf &Quaternionf::Set(float _x, float _y, float _z, float _w)
{
	x = _x;
	y = _y;
	z = _z;
	w = _w;

	return (*this);
}

inline Quaternionf &Quaternionf::Set(float vector[4])
{
	x = vector[0];
	y = vector[1];
	z = vector[2];
	w = vector[3];

	return (*this);
}

inline Quaternionf &Quaternionf::Set(const Vec4f &vec4)
{
	x = vec4.x;
	y = vec4.y;
	z = vec4.z;
	w = vec4.w;
}

float Quaternionf::GetAngle() const
{
	double angle = 2.0 * acos(w);

	return static_cast<float>(angle <= RT_PI? angle : RT_2PI - angle);
}

Mat4f Quaternionf::Matrix4()
{
	float a00 = 2.f * x * x;
	float a11 = 2.f * y * y;
	float a22 = 2.f * z * z;

	float a01 = 2.f * x * y;
	float a02 = 2.f * x * z;
	float a03 = 2.f * x * w;

	float a12 = 2.f * y * z;
	float a13 = 2.f * y * w;

	float a23 = 2.f * z * w;

	return Mat4f(1.f - a11 - a22, a01 + a23      , a02 - a13      , 0.f,
				  a01 - a23      , 1.f - a22 - a00, a12 + a03      , 0.f,
				  a02 + a13      , a12 - a03      , 1.f - a11 - a00, 0.f,
				  0.f            , 0.f            , 0.f            , 1.0f);
}

Quaternionf &Quaternionf::operator=(const Quaternionf &other)
{
	if (this != &other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
		w = other.w;
	}

	return (*this);
}

Quaternionf &Quaternionf::operator+=(const Quaternionf &other)
{
	x += other.x;
	y += other.y;
	z += other.z;
	w += other.w;

	return (*this);
}

Quaternionf Quaternionf::operator+(const Quaternionf &other) const
{
	return Quaternionf(x + other.x,
					    y + other.y,
					    z + other.z,
					    w + other.w);
}

Quaternionf &Quaternionf::operator-=(const Quaternionf &other)
{
	x -= other.x;
	y -= other.y;
	z -= other.z;
	w -= other.w;

	return (*this);
}

Quaternionf Quaternionf::operator-(const Quaternionf &other) const
{
	return Quaternionf(x - other.x,
					    y - other.y,
					    z - other.z,
					    w - other.w);
}

Quaternionf &Quaternionf::operator*=(const Quaternionf &other)
{
	x = w * other.x + x * other.w + y * other.z - z * other.y;
	y = w * other.y - x * other.z + y * other.w + z * other.x;
	z = w * other.z + x * other.y - y * other.x + z * other.w;
	w = w * other.w - x * other.x - y * other.y - z * other.z;

	return (*this);
}

Quaternionf &Quaternionf::operator*=(float scalar)
{
	x *= scalar;
	y *= scalar;
	z *= scalar;
	w *= scalar;

	return (*this);
}

Quaternionf Quaternionf::operator*(const Quaternionf &other) const
{
	return Quaternionf(w * other.x + x * other.w + y * other.z - z * other.y,
					    w * other.y - x * other.z + y * other.w + z * other.x,
					    w * other.z + x * other.y - y * other.x + z * other.w,
					    w * other.w - x * other.x - y * other.y - z * other.z);
}

Quaternionf Quaternionf::operator*(float scalar) const
{
	return Quaternionf(x*scalar, y*scalar,
					    z*scalar, w*scalar);
}

Quaternionf &Quaternionf::operator/=(const Quaternionf &other)
{
	return (*this *= Qtf::Inverse(other));
}

Quaternionf &Quaternionf::operator/=(float scalar)
{
	x /= scalar;
	y /= scalar;
	z /= scalar;
	w /= scalar;

	return (*this);
}

Quaternionf Quaternionf::operator/(const Quaternionf &other) const
{
	return (*this * Qtf::Inverse(other));
}

Quaternionf Quaternionf::operator/(float scalar) const
{
	return Quaternionf(x/scalar, y/scalar,
					    z/scalar, w/scalar);
}

Quaternionf Quaternionf::operator-() const
{
	return Quaternionf(-x, -y, -z, -w);
}

Quaternionf &Quaternionf::Normalize()
{
	float nor = 1.0f / Size();

	x *= nor;
	y *= nor;
	z *= nor;
	w *= nor;

	return (*this);
}

//Quaternionf Quaternionf::Inverse()
//{
//	return Qtf::Conjugate(*this) / SizeSQR();
//}

Quaternionf &Quaternionf::Slerp(const Quaternionf &target, float alpha)
{
	float cosO = Qtf::asVec4(*this).Dot(Qtf::asVec4(target));
	float absCosO = abs(cosO);
	float scale1, scale2;

	if (1.f - absCosO > 1E-6f)
	{
		float sinSQR = 1.f - absCosO * absCosO;
		float sinO = (1.f / sqrtf(sinSQR));
		float o = atan2(sinSQR * sinO, absCosO);

		scale1 = sinf((1.f - alpha) * o) * sinO;
		scale2 = sinf(alpha * o) * sinO;
	}
	else
	{
		scale1 = 1.f - alpha;
		scale2 = alpha;
	}

	scale2 = cosO >= 0.f ? scale1 : -scale1;
	x = scale1 * x + scale2 * target.x;
	y = scale1 * y + scale2 * target.y;
	z = scale1 * z + scale2 * target.z;
	w = scale1 * w + scale2 * target.w;

	return (*this);
}

Quaternionf &Quaternionf::Nlerp(const Quaternionf &q, float factor)
{
	float cosO = Qtf::asVec4(*this).Dot(Qtf::asVec4(q));
	float scale1 = 1.f - factor;
	float scale2 = (cosO >= 0.f) ? factor : -factor;

	x = scale1 * x + scale2 * q.x;
	y = scale1 * y + scale2 * q.y;
	z = scale1 * z + scale2 * q.z;
	w = scale1 * w + scale2 * q.w;

	return Normalize();
}

Quaternionf &Quaternionf::Rotate(float angleX, float angleY, float angleZ)
{
	Vec3f theta = Vec3f(angleX * 0.5f, 
						  angleY * 0.5f,
						  angleZ * 0.5f);

	float lenght = theta.SizeSQR();
	float qx, qy, qz, qw, s;

	if (lenght * lenght / 24.f < 1E-8f)
	{
		qw = 1.f - lenght / 2.f;
		s = 1.f - lenght / 6.f;
	}
	else
	{
		float tm = theta.Size();
		qw = cosf(tm);
		s = sinf(tm) / tm;
	}

	qx = theta.x * s;
	qy = theta.y * s;
	qz = theta.z * s;

	Set(w * qx + x * qw + y * qz - z * qy,
		w * qy - x * qz + y * qw + z * qx,
		w * qz + x * qy - y * qx + z * qw,
		w * qw - x * qx - y * qy - z * qz);

	return (*this);
}

Quaternionf &Quaternionf::RotateX(float angle)
{
	float Cos = cos(angle * 0.5f);
	float Sin = sin(angle * 0.5f);

	Set(w * Sin + x * Cos,
		y * Cos + z * Sin,
		z * Cos - y * Sin,
		w * Cos - x * Sin);

	return (*this);
}

Quaternionf &Quaternionf::RotateY(float angle)
{
	float Cos = cos(angle * 0.5f);
	float Sin = sin(angle * 0.5f);

	Set(x * Cos - z * Sin,
		w * Sin + y * Cos,
		x * Sin + z * Cos,
		w * Cos - y * Sin);

	return (*this);
}

Quaternionf &Quaternionf::RotateZ(float angle)
{
	float Cos = cos(angle * 0.5f);
	float Sin = sin(angle * 0.5f);

	Set(x * Cos + y * Sin,
		y * Cos - x * Sin,
		w * Sin + z * Cos,
		w * Cos - z * Sin);

	return (*this);
}

inline float Quaternionf::Size() const
{
	float length = x*x + y*y + z*z + w*w;
	return sqrtf(length);
}

inline float Quaternionf::SizeSQR() const
{
	return x*x + y*y + z*z + w*w;
}


