#ifndef __RT_MAT4_H__
#define __RT_MAT4_H__

#include <ostream>
#include <iostream>
#include "Util.h"
#include "RTvec3.h"
#include "RTvec4.h"

namespace RT
{
	class Mat3f;

	class Mat4f
	{
	private:
		float matrix[4][4];

	public:
		explicit Mat4f(float m0x0 = 1, float m0x1 = 0, float m0x2 = 0, float m0x3 = 0,
					   float m1x0 = 0, float m1x1 = 1, float m1x2 = 0, float m1x3 = 0,
					   float m2x0 = 0, float m2x1 = 0, float m2x2 = 1, float m2x3 = 0,
					   float m3x0 = 0, float m3x1 = 0, float m3x2 = 0, float m3x3 = 1);

		Mat4f(const Mat4f &other);
		explicit Mat4f(float matrix4x4[SIZE_4][SIZE_4]);
		explicit Mat4f(Mat3f &matrix3x3);
		~Mat4f();

		Mat4f &Set(float m0x0, float m0x1, float m0x2, float m0x3,
				   float m1x0, float m1x1, float m1x2, float m1x3,
				   float m2x0, float m2x1, float m2x2, float m2x3,
				   float m3x0, float m3x1, float m3x2, float m3x3);

		Mat4f &Set(const Mat4f &other);
		Mat4f &Set(float matrix4x4[SIZE_4][SIZE_4]);
		Mat4f &Set(Mat3f &matrix3x3);
		Mat4f &Set(float vector[SIZE_4]);
		Mat4f &Set(const Vec4f &vector4);
		Mat4f &Set(int i, int j, float value);

		inline float *GetMatrix() { return matrix[0]; }
		inline float Getixj(int i, int j) const { return matrix[i][j]; }
		inline int GetWidth() const { return SIZE_4; }
		inline int GetHeight() const { return SIZE_4; }
		inline int Length() const { return SIZE_4 * SIZE_4; }

		Mat4f &Transpose();
		Mat4f &Inverse();
		Mat4f &Identity();
		Mat4f &Opposed();

		Mat4f &Null();

		float Determinant() const;

		inline float *operator[](const int index) { return matrix[index]; }
		inline float *operator+(const int index) { return matrix[index]; }

		Mat4f &operator=(const Mat4f &other);
		Mat4f &operator=(Mat3f &matrix3x3);
		bool operator==(const Mat4f &other) const;
		bool operator!=(const Mat4f &other) const;

		Mat4f operator+(const Mat4f &other) const;
		Mat4f operator+(float matrix[SIZE_4][SIZE_4]) const;
		Mat4f operator+(float scalar) const;
		Mat4f &operator+=(const Mat4f &other);
		Mat4f &operator+=(float matrix[SIZE_4][SIZE_4]);
		Mat4f &operator+=(float scalar);

		Mat4f operator-(const Mat4f &other) const;
		Mat4f operator-(float matrix[SIZE_4][SIZE_4]) const;
		Mat4f operator-(float scalar) const;
		Mat4f &operator-=(const Mat4f &other);
		Mat4f &operator-=(float matrix[SIZE_4][SIZE_4]);
		Mat4f &operator-=(float scalar);

		Mat4f operator*(const Mat4f &other) const;
		Mat4f operator*(Mat3f &matrix) const;
		Mat4f operator*(float matrix[SIZE_4][SIZE_4]) const;
		Vec4f operator*(const Vec4f &vector) const;
		Mat4f operator*(float scalar) const;
		Mat4f &operator*=(const Mat4f &other);
		Mat4f &operator*=(Mat3f &matrix);
		Mat4f &operator*=(float matrix[SIZE_4][SIZE_4]);
		Mat4f &operator*=(float scalar);

		Mat4f operator/(float scalar) const;
		Mat4f &operator/=(float scalar);

		//Transformações
		Mat4f &Scale(const Vec3f &vector);
		Mat4f &Scale(float x, float y, float z);
		Mat4f &Scale(float size);

		Mat4f &RotateX(float angle);
		Mat4f &RotateY(float angle);
		Mat4f &RotateZ(float angle);
		Mat4f &RotateXYZ(float angleX, float angleY, float angleZ);

		Mat4f &Translate(const Vec3f &distance);
		Mat4f &Translate(float x, float y, float z);

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		Mat4f &Ortho(const float left, const float right, const float bottom,
			const float top, const float near, const float far);

		Mat4f &Perspective(const Vec4f &parameters);
		Mat4f &Perspective(const float fovy, const float aspect,
			const float near, const float far);

		Mat4f &LookAt(const Vec3f &eye, const Vec3f &center, const Vec3f &up);
		Mat4f &LookAt(const float xEye, const float yEye, const float zEye,
					  const float xCenter, const float yCenter, const float zCenter,
					  const float xUp, const float yUp, const float zUp);

		Vec4f Project(const Vec3f *vector, const int *viewport);
		Vec4f Project(const float x, const float y, const float z, const int *viewport);
		Vec3f Project3f(const Vec3f *vector, const int *viewport);
		Vec3f Project3f(const float x, const float y, const float z, const int *viewport);

		friend std::ostream &operator<<(std::ostream &out, Mat4f &mat)
		{
			for (int i = 0; i < SIZE_4; i++)
			{
				for (int j = 0; j < SIZE_4; j++)
					out << mat.matrix[i][j] << " ";
				out << std::endl;
			}

			std::cout << std::endl;
			return out;
		}
	};

	Mat4f operator*(float scalar, const Mat4f &matrix);
	Mat4f operator+(float scalar, const Mat4f &matrix);
	Mat4f operator-(float scalar, const Mat4f &matrix);

	namespace mt4
	{
		inline Mat4f Transpose(const Mat4f &mat)
		{
			return Mat4f(mat).Transpose();
		}

		inline Mat4f Inverse(const Mat4f &mat)
		{
			return Mat4f(mat).Inverse();
		}

		inline Mat4f Opposed(const Mat4f &mat)
		{
			return Mat4f(mat).Opposed();
		}

		inline Mat4f Ortho(const float left, const float right, const float bottom,
						   const float top, const float near, const float far)
		{
			return Mat4f().Ortho(left, right, bottom, top, near, far);
		}

		inline  Mat4f Perspective(const Vec4f &parameters)
		{
			return Mat4f().Perspective(parameters);
		}

		inline Mat4f Perspective(const float fovy, const float aspect,
								 const float near, const float far)
		{
			return Mat4f().Perspective(fovy, aspect, near, far);
		}

		inline Mat4f LookAt(const Vec3f &eye, const Vec3f &center, const Vec3f &up)
		{
			return Mat4f().LookAt(eye, center, up);
		}

		inline Mat4f LookAt(const float xEye, const float yEye, const float zEye,
							const float xCenter, const float yCenter, const float zCenter,
							const float xUp, const float yUp, const float zUp)
		{
			return Mat4f().LookAt(xEye, yEye, zEye,
				xCenter, yCenter, zCenter,
				xUp, yUp, zUp);
		}


		inline Mat4f AffTranslation(const Vec3f &distance)
		{
			return Mat4f(1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				distance.x, distance.y, distance.z, 1);
		}

		inline Mat4f AffTranslation(float x, float y, float z)
		{
			return Mat4f(1, 0, 0, 0,
						 0, 1, 0, 0,
						 0, 0, 1, 0,
						 x, y, z, 1);
		}

		inline Mat4f AffRotation(float angleX, float angleY, float angleZ)
		{
			return Mat4f().RotateXYZ(angleX, angleY, angleZ);
		}

		inline Mat4f AffRotation(const Vec3f &angles)
		{
			return Mat4f().RotateXYZ(angles.x, angles.y, angles.z);
		}

		inline Mat4f AffRotationX(float angle)
		{
			float Cos = cos(angle);
			float Sin = sin(angle);

			return Mat4f(1, 0, 0, 0,
						 0, Cos, Sin, 0,
						 0, -Sin, Cos, 0,
						 0, 0, 0, 1);
		}

		inline Mat4f AffRotationY(float angle)
		{
			float Cos = cos(angle);
			float Sin = sin(angle);

			return Mat4f(Cos, 0, -Sin, 0,
						 0, 1, 0, 0,
						 Sin, 0, Cos, 0,
						 0, 0, 0, 1);
		}

		inline Mat4f AffRotationZ(float angle)
		{
			float Cos = cos(angle);
			float Sin = sin(angle);

			return Mat4f(Cos, Sin, 0, 0,
						-Sin, Cos, 0, 0,
						0, 0, 1, 0,
						0, 0, 0, 1);
		}

		inline Mat4f AffScaling(const Vec3f &scales)
		{
			return Mat4f(scales.x, 0, 0, 0,
				0, scales.y, 0, 0,
				0, 0, scales.z, 0,
				0, 0, 0, 1);
		}

		inline Mat4f AffScaling(float scaleX, float scaleY, float scaleZ)
		{
			return Mat4f(scaleX, 0, 0, 0,
				0, scaleY, 0, 0,
				0, 0, scaleZ, 0,
				0, 0, 0, 1);
		}

		inline Mat4f AffScaling(float scale)
		{
			return AffScaling(scale, scale, scale);
		}

		inline Vec3f Transform(const Mat4f &matrix, const Vec3f &vector, float w)
		{
			return (matrix * Vec4f(vector, w)).XYZ();
		}

		inline Vec4f Transform(const Mat4f &matrix, const Vec4f &vector4)
		{
			return (matrix * vector4);
		}
	}
}
#endif //__RT_MAT4_H__
