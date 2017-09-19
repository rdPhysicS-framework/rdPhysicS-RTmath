#ifndef __RT_MAT3_H__
#define __RT_MAT3_H__

#include "Util.h"
#include <ostream>

#include "RTvec2.h"
#include "RTvec3.h"

namespace RT
{
	class Mat4f;

	class Mat3f
	{
	private:
		float matrix[3][3];

	public:
		explicit Mat3f(float m0x0 = 1, float m0x1 = 0, float m0x2 = 0,
					   float m1x0 = 0, float m1x1 = 1, float m1x2 = 0,
					   float m2x0 = 0, float m2x1 = 0, float m2x2 = 1);
		Mat3f(const Mat3f &matrix);
		explicit Mat3f(float matrix[SIZE_3][SIZE_3]);
		~Mat3f();

		Mat3f &Set(float m0x0, float m0x1, float m0x2,
				   float m1x0, float m1x1, float m1x2,
				   float m2x0, float m2x1, float m2x2);
		Mat3f &Set(const Mat3f &other);
		Mat3f &Set(float _matrix[SIZE_3][SIZE_3]);
		Mat3f &Set(const Vec3f &vector);
		Mat3f &Set(const Vec2f &vector);

		/*void SetPtr(float ptr[9]);

		RTmatrix *GetMatrix() const;*/
		inline float *GetMatrix() { return matrix[0]; }
		float Getixj(int i, int j) const;
		int GetWidth() const;
		int GetHeight() const;
		int Length() const;

		Mat3f &Transpose();
		Mat3f &Inverse();
		Mat3f &Identity();
		Mat3f &Opposed();
		Mat3f &Null();
		float Determinant() const;

		Mat3f &operator=(const Mat3f &other);
		inline float *operator[](const int index) { return matrix[index]; }
		inline float *operator+(const int index) { return matrix[index]; }

		bool operator==(const Mat3f &other) const;
		bool operator!=(const Mat3f &other) const;

		Mat3f &operator+=(const Mat3f &other);
		Mat3f &operator+=(float matrix[SIZE_3][SIZE_3]);
		Mat3f &operator+=(float scalar);
		Mat3f operator+(const Mat3f &other) const;
		Mat3f operator+(float matrix[SIZE_3][SIZE_3]) const;
		Mat3f operator+(float scalar) const;

		Mat3f &operator-=(const Mat3f &other);
		Mat3f &operator-=(float matrix[SIZE_3][SIZE_3]);
		Mat3f &operator-=(float scalar);
		Mat3f operator-(const Mat3f &other) const;
		Mat3f operator-(float matrix[SIZE_3][SIZE_3]) const;
		Mat3f operator-(float scalar) const;

		Mat3f &operator*=(const Mat3f &other);
		Mat3f &operator*=(Mat4f &matrix);
		Mat3f &operator*=(float matrix[SIZE_3][SIZE_3]);
		Mat3f operator*(const Mat3f &other) const;
		Mat3f operator*(Mat4f &matrix) const;
		Mat3f operator*(float matrix[SIZE_3][SIZE_3]) const;
		Vec3f operator*(const Vec3f &vector) const;
		Mat3f &operator*=(float scalar);
		Mat3f operator*(float scalar) const;

		Mat3f &operator/=(float scalar);
		Mat3f operator/(float scalar) const;

		Mat3f &Translate(float _x, float _y);
		Mat3f &Translate(const Vec2f &vector);

		Mat3f &Scale(float scale);
		Mat3f &Scale(float scalex, float scaley);
		Mat3f &Scale(const Vec2f &scales);

		Mat3f &Rotate(float angle);
	};

	Mat3f operator*(float scalar, const Mat3f &matrix);
	Mat3f operator+(float scalar, const Mat3f &matrix);
	Mat3f operator-(float scalar, const Mat3f &matrix);
	std::ostream &operator<<(std::ostream &output, Mat3f &mat);

	namespace mt3
	{
		inline Mat3f AffTranslation(const Vec2f &distance)
		{
			return Mat3f(1.0f, 0.0f, distance.x,
						 0.0f, 1.0f, distance.y,
						 0.0f, 0.0f, 1.0f);
		}

		inline Mat3f AffTranslation(float x, float y)
		{
			return Mat3f(1.0f, 0.0f, x,
						 0.0f, 1.0f, y,
						 0.0f, 0.0f, 1.f);
		}

		inline Mat3f AffRotation(float angle)
		{
			float Cos = cos(angle);
			float Sin = sin(angle);

			return Mat3f(Cos, -Sin, 0.0f,
						 Sin, Cos,  0.0f,
						 0.0f, 0.0f,  1.0f);
		}

		inline Mat3f AffScaling(const Vec2f &scales)
		{
			return Mat3f(scales.x, 0.0f,     0.0f,
						 0.0f,     scales.y, 0.0f,
						 0.0f,     0.0f,     1.0f);
		}

		inline Mat3f AffScaling(float scaleX, float scaleY)
		{
			return Mat3f(scaleX, 0.f,    0.f,
				         0.f,    scaleY, 0.f,
				         0.f,    0.f,    1.f);
		}

		inline Mat3f AffScaling(float scale)
		{
			return AffScaling(scale, scale);
		}

		inline Vec2f Transform(const Mat3f &matrix, const Vec2f &vector2, float z = 1.0f)
		{
			return (matrix * Vec3f(vector2, z)).XY();
		}

		inline Vec3f Transform(const Mat3f &matrix, const Vec3f &vector3)
		{
			return matrix * vector3;
		}

		inline Mat3f Transpose(const Mat3f &mat)
		{
			return Mat3f(mat).Transpose();
		}

		inline Mat3f Inverse(const Mat3f &mat)
		{
			return Mat3f(mat).Inverse();
		}

		inline Mat3f Identity(const Mat3f &mat)
		{
			return Mat3f(mat).Identity();
		}

		inline Mat3f Opposed(const Mat3f &mat)
		{
			return Mat3f(mat).Opposed();
		}
	}
}
#endif //__RT_MAT3_H__