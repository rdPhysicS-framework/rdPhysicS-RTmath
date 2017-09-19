#ifndef __RT_QUATERNION_H__ 
#define __RT_QUATERNION_H__ 

#include "Util.h"
#include <sstream>
#include <ostream>

#include "RTvec3.h"
#include "RTvec4.h"

namespace RT
{
	class Mat4f;

	class Quaternionf
	{
	public:
		float x;
		float y;
		float z;
		float w;

	public:
		explicit Quaternionf(float _x = 0, float _y = 0, float _z = 0, float _w = 1);
		explicit Quaternionf(float vector[SIZE_4]);
		explicit Quaternionf(const Vec4f &vec4);
		explicit Quaternionf(float angle, const Vec3f &axis);
		explicit Quaternionf(float angle1, const Vec3f &axis1,
			float angle2, const Vec3f &axis2,
			float angle3, const Vec3f &axis3);
		Quaternionf(const Quaternionf &other);
		~Quaternionf();

		Quaternionf &RotationAxis(float angle, const Vec3f &axis);
		Quaternionf &RotationAxis(float angle, float _x,
			float _y, float _z);
		Quaternionf &RotationAxis(float angle1, const Vec3f &axis1,
			float angle2, const Vec3f &axis2,
			float angle3, const Vec3f &axis3);

		inline Quaternionf &Set(float _x, float _y, float _z, float _w);
		inline Quaternionf &Set(float vector[4]);
		inline Quaternionf &Set(const Vec4f &vec4);

		float GetAngle() const;

		Mat4f Matrix4();

		inline float operator[](int index) const
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;
			if (index == 3) return w;

			std::stringstream s;
			s << "Index out of bounds " << index << ", size 3.\n";

			throw std::out_of_range(s.str());
		}

		inline float &operator[](int index)
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;
			if (index == 3) return w;

			std::stringstream s;
			s << "Index out of bounds " << index << ", size 3.\n";

			throw std::out_of_range(s.str());
		}

		inline float operator+(int index) const
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;
			if (index == 3) return w;

			std::stringstream s;
			s << "Index out of bounds " << index << ", size 3.\n";

			throw std::out_of_range(s.str());
		}

		inline float &operator+(int index)
		{
			if (index == 0) return x;
			if (index == 1) return y;
			if (index == 2) return z;
			if (index == 3) return w;

			std::stringstream s;
			s << "Index out of bounds " << index << ", size 3.\n";

			throw std::out_of_range(s.str());
		}

		Quaternionf &operator=(const Quaternionf &other);
		inline bool operator==(const Quaternionf &other) const
		{
			return (x == other.x) &&
				(y == other.y) &&
				(z == other.z) &&
				(w == other.w);
		}
		inline bool operator!=(const Quaternionf &other) const
		{
			return (x != other.x) ||
				(y != other.y) ||
				(z != other.z) ||
				(w != other.w);
		}

		Quaternionf &operator+=(const Quaternionf &other);
		Quaternionf operator+(const Quaternionf &other) const;
		Quaternionf &operator-=(const Quaternionf &other);
		Quaternionf operator-(const Quaternionf &other) const;
		Quaternionf &operator*=(const Quaternionf &other);
		Quaternionf &operator*=(float scalar);
		Quaternionf operator*(const Quaternionf &other) const;
		Quaternionf operator*(float scalar) const;
		Quaternionf &operator/=(const Quaternionf &other);
		Quaternionf &operator/=(float scalar);
		Quaternionf operator/(const Quaternionf &other) const;
		Quaternionf operator/(float scalar) const;
		Quaternionf operator-() const;

		Quaternionf &Normalize();
		//Quaternionf Inverse();
		Quaternionf &Slerp(const Quaternionf &target, float alpha);
		Quaternionf &Nlerp(const Quaternionf &q, float factor);

		Quaternionf &Rotate(float angleX, float angleY, float angleZ);
		Quaternionf &RotateX(float angle);
		Quaternionf &RotateY(float angle);
		Quaternionf &RotateZ(float angle);

		inline float Size() const;
		inline float SizeSQR() const;

		friend std::ostream &operator<<(std::ostream &out, const Quaternionf &q)
		{
			out << "(" << q.x << ", " << q.y << ", "
				<< q.z << ", " << q.w << ")\n";
			return out;
		}
	};

	inline Quaternionf operator*(float scalar, const Quaternionf &q)
	{
		return (q * scalar);
	}


	namespace Qtf
	{
		inline Vec3f asVec3(const Quaternionf &q)
		{
			return Vec3f(q.x, q.y, q.z);
		}

		inline Vec4f asVec4(const Quaternionf &q)
		{
			return Vec4f(q.x, q.y, q.z, q.w);
		}

		inline Quaternionf Conjugate(const Quaternionf &q)
		{
			return Quaternionf(-q.x, -q.y, -q.z, q.w);
		}

		inline Quaternionf Inverse(const Quaternionf &q)
		{
			return Qtf::Conjugate(q) / q.SizeSQR();
		}

		inline Quaternionf Slerp(float alpha, const Quaternionf &from, Quaternionf &to)
		{
			Quaternionf aux(to);

			double cosO = Qtf::asVec4(from).Dot(Qtf::asVec4(to));
			double epsilo = 1E-6f;

			if (cosO < epsilo)
			{
				cosO *= -1;
				aux = -to;
			}

			double o;
			double sinO;
			double scale1, scale2;

			if (1.0 - cosO > epsilo)
			{
				o = acos(cosO);
				sinO = sin(o);

				scale1 = sin((1.0 - alpha) * o) / sinO;
				scale2 = sin(alpha * o) / sinO;
			}
			else
			{
				scale1 = 1.0 - alpha;
				scale2 = alpha;
			}

			return Quaternionf(((from * scale1) + (aux * scale2)));
		}

		inline Quaternionf Nlerp(float factor, const Quaternionf &from, Quaternionf &to)
		{
			float cosO = Qtf::asVec4(from).Dot(Qtf::asVec4(to));
			float scale1 = 1.f - factor;
			float scale2 = (cosO >= 0.f) ? factor : -factor;

			Quaternionf aux = Quaternionf(scale1 * from.x + scale2 * to.x,
											scale1 * from.y + scale2 * to.y,
											scale1 * from.z + scale2 * to.z,
											scale1 * from.w + scale2 * to.w);

			aux.Normalize();

			return aux;
		}

		inline Quaternionf Normlize(const Quaternionf &q)
		{
			return Quaternionf(q).Normalize();
		}

		inline Quaternionf Rotate(const Quaternionf &q, float angleX, float angleY, float angleZ)
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

			return Quaternionf(q.w * qx + q.x * qw + q.y * qz - q.z * qy,
							    q.w * qy - q.x * qz + q.y * qw + q.z * qx,
							    q.w * qz + q.x * qy - q.y * qx + q.z * qw,
							    q.w * qw - q.x * qx - q.y * qy - q.z * qz);
		}

		inline Quaternionf RotateAxis(const Quaternionf &q, float angle, float axisX, float axisy, float axisZ)
		{
			double a = angle / 2.0;
			double sina = sin(a);

			Vec3f aux = vc3::Normalize(Vec3f(axisX, axisy, axisZ)) * sina;
			double w = cos(a);

			return Quaternionf(q.w * aux.x + q.x * w + q.y * aux.z - q.z * aux.y,
								q.w * aux.y - q.x * aux.z + q.y * w + q.z * aux.x,
								q.w * aux.z + q.x * aux.y - q.y * aux.x + q.z * w,
								q.w * w - q.x * aux.x - q.y * aux.y - q.z * aux.z);
		}

		inline Quaternionf RotateAxis(const Quaternionf &q, float angle, const Vec3f &axis)
		{
			double a = angle / 2.0;
			double sina = sin(a);

			Vec3f aux = vc3::Normalize(axis) * sina;
			double w = cos(a);

			return Quaternionf(q.w * aux.x + q.x * w + q.y * aux.z - q.z * aux.y,
								q.w * aux.y - q.x * aux.z + q.y * w + q.z * aux.x,
								q.w * aux.z + q.x * aux.y - q.y * aux.x + q.z * w,
								q.w * w - q.x * aux.x - q.y * aux.y - q.z * aux.z);
		}
	}
}
#endif //__RT_QUATERNION_H__ 

