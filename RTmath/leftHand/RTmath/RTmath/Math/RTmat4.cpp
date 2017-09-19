#include "RTmat4.h"
#include "RTmat3.h"
#include "RTmathOp.h"
#include<cmath>
#include<iostream>

using namespace RT;

Mat4f::Mat4f(float m0x0, float m0x1, float m0x2, float m0x3,
			   float m1x0, float m1x1, float m1x2, float m1x3,
			   float m2x0, float m2x1, float m2x2, float m2x3, 
			   float m3x0, float m3x1, float m3x2, float m3x3) 
{
	Set(m0x0, m0x1, m0x2, m0x3,
		m1x0, m1x1, m1x2, m1x3,
		m2x0, m2x1, m2x2, m2x3,
		m3x0, m3x1, m3x2, m3x3 );
}

Mat4f::Mat4f(const Mat4f &other) 
{
	memcpy(matrix, (void*)other.matrix, sizeof(matrix));
}

Mat4f::Mat4f(float matrix4x4[SIZE_4][SIZE_4])
{
	memcpy(matrix, (void*)matrix4x4, sizeof(matrix));
}

Mat4f::Mat4f(Mat3f &matrix3x3)
{
	Set(matrix3x3);
}

Mat4f::~Mat4f()
{}

Mat4f &Mat4f::Set(float m0x0, float m0x1, float m0x2, float m0x3, 
				    float m1x0, float m1x1, float m1x2, float m1x3, 
					float m2x0, float m2x1, float m2x2, float m2x3,
					float m3x0, float m3x1, float m3x2, float m3x3)
{
	matrix[0][0] = m0x0; matrix[0][1] = m0x1; matrix[0][2] = m0x2; matrix[0][3] = m0x3;
	matrix[1][0] = m1x0; matrix[1][1] = m1x1; matrix[1][2] = m1x2; matrix[1][3] = m1x3;
	matrix[2][0] = m2x0; matrix[2][1] = m2x1; matrix[2][2] = m2x2; matrix[2][3] = m2x3;
	matrix[3][0] = m3x0; matrix[3][1] = m3x1; matrix[3][2] = m3x2; matrix[3][3] = m3x3;

	return (*this);
}

Mat4f &Mat4f::Set(const Mat4f &other)
{
	memcpy(matrix, (void*)other.matrix, sizeof(matrix));

	return (*this);
}

Mat4f &Mat4f::Set(float matrix4x4[SIZE_4][SIZE_4])
{
	memcpy(matrix, (void*)matrix4x4, sizeof(matrix));

	return (*this);
}

Mat4f &Mat4f::Set(Mat3f &matrix3x3)
{
	matrix[0][0] = matrix3x3[0][0]; 
	matrix[0][1] = matrix3x3[0][1]; 
	matrix[0][2] = matrix3x3[0][2];
	matrix[0][3] = 0;

	matrix[1][0] = matrix3x3[1][0];
	matrix[1][1] = matrix3x3[1][1]; 
	matrix[1][2] = matrix3x3[1][2];
	matrix[1][3] = 0;

	matrix[2][0] = matrix3x3[2][0]; 
	matrix[2][1] = matrix3x3[2][1]; 
	matrix[2][2] = matrix3x3[2][2];
	matrix[2][3] = 0;

	matrix[3][0] = 0;
	matrix[3][1] = 0;
	matrix[3][2] = 0;
	matrix[3][3] = 1;

	return (*this);
}

Mat4f &Mat4f::Set(float vector[SIZE_4])
{
	matrix[3][0] = vector[0];
	matrix[3][1] = vector[1];
	matrix[3][2] = vector[2];
	matrix[3][3] = vector[3];


	return (*this);
}

Mat4f &Mat4f::Set(const Vec4f &vector4)
{
	matrix[3][0] = vector4.x;
	matrix[3][1] = vector4.y;
	matrix[3][2] = vector4.z;
	matrix[3][3] = vector4.w;

	return (*this);
}

Mat4f &Mat4f::Set(int i, int j, float value)
{
	matrix[i][j] = value;

	return (*this);
}

Mat4f &Mat4f::Transpose()
{
	Mat4f aux(*this);

	matrix[0][0] = aux.matrix[0][0];
	matrix[0][1] = aux.matrix[1][0];
	matrix[0][2] = aux.matrix[2][0];
	matrix[0][3] = aux.matrix[3][0];

	matrix[1][0] = aux.matrix[0][1];
	matrix[1][1] = aux.matrix[1][1];
	matrix[1][2] = aux.matrix[2][1];
	matrix[1][3] = aux.matrix[3][1];

	matrix[2][0] = aux.matrix[0][2];
	matrix[2][1] = aux.matrix[1][2];
	matrix[2][2] = aux.matrix[2][2];
	matrix[2][3] = aux.matrix[3][2];

	matrix[3][0] = aux.matrix[0][3];
	matrix[3][1] = aux.matrix[1][3];
	matrix[3][2] = aux.matrix[2][3];
	matrix[3][3] = aux.matrix[3][3];

	return (*this);
}

Mat4f &Mat4f::Inverse()
{
	float a = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	float b = matrix[0][0] * matrix[1][2] - matrix[0][2] * matrix[1][0];
	float c = matrix[0][0] * matrix[1][3] - matrix[0][3] * matrix[1][0];
	float d = matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1];
	float e = matrix[0][1] * matrix[1][3] - matrix[0][3] * matrix[1][1];
	float f = matrix[0][2] * matrix[1][3] - matrix[0][3] * matrix[1][2];
	float g = matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0];
	float h = matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0];
	float i = matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0];
	float j = matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1];
	float k = matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1];
	float l = matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2];
	float det = a * l - b * k + c * j + d * i - e * h + f * g;
	det = 1.0f / det;

	return Set((matrix[1][1] * l - matrix[1][2] * k + matrix[1][3] * j) * det,
			  (-matrix[0][1] * l + matrix[0][2] * k - matrix[0][3] * j) * det,
			   (matrix[3][1] * f - matrix[3][2] * e + matrix[3][3] * d) * det,
			  (-matrix[2][1] * f + matrix[2][2] * e - matrix[2][3] * d) * det,
			  (-matrix[1][0] * l + matrix[1][2] * i - matrix[1][3] * h) * det,
			   (matrix[0][0] * l - matrix[0][2] * i + matrix[0][3] * h) * det,
			  (-matrix[3][0] * f + matrix[3][2] * c - matrix[3][3] * b) * det,
			   (matrix[2][0] * f - matrix[2][2] * c + matrix[2][3] * b) * det,
			   (matrix[1][0] * k - matrix[1][1] * i + matrix[1][3] * g) * det,
			  (-matrix[0][0] * k + matrix[0][1] * i - matrix[0][3] * g) * det,
			   (matrix[3][0] * e - matrix[3][1] * c + matrix[3][3] * a) * det,
			  (-matrix[2][0] * e + matrix[2][1] * c - matrix[2][3] * a) * det,
			  (-matrix[1][0] * j + matrix[1][1] * h - matrix[1][2] * g) * det,
			   (matrix[0][0] * j - matrix[0][1] * h + matrix[0][2] * g) * det,
			  (-matrix[3][0] * d + matrix[3][1] * b - matrix[3][2] * a) * det,
			   (matrix[2][0] * d - matrix[2][1] * b + matrix[2][2] * a) * det);

	return (*this);
}

Mat4f &Mat4f::Identity()
{
	return Set(1.0f, 0.0f, 0.0f, 0.0f,
			   0.0f, 1.0f, 0.0f, 0.0f,
			   0.0f, 0.0f, 1.0f, 0.0f,
			   0.0f, 0.0f, 0.0f, 1.0f);
}

Mat4f &Mat4f::Opposed()
{
	return (*this *= -1);
}

Mat4f &Mat4f::Null()
{
	matrix[0][0] = matrix[0][1] = matrix[0][2] =  matrix[0][3] = 0;
	matrix[1][0] = matrix[1][1] = matrix[1][2] =  matrix[1][3] = 0;
	matrix[2][0] = matrix[2][1] = matrix[2][2] =  matrix[2][3] = 0;
	matrix[3][0] = matrix[3][1] = matrix[3][2] =  matrix[3][3] = 0;

	return (*this);
}

float Mat4f::Determinant() const
{
	return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) *
		   (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) +

		   (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) *
		   (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) +

		   (matrix[0][0] * matrix[1][3] - matrix[0][3] * matrix[1][0]) *
		   (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1]) +

		   (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) *
		   (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +

		   (matrix[0][3] * matrix[1][1] - matrix[0][1] * matrix[1][3]) *
		   (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0]) +

		   (matrix[0][2] * matrix[1][3] - matrix[0][3] * matrix[1][2]) *
		   (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0]);
}

Mat4f &Mat4f::operator=(const Mat4f &other)
{
	if (this != &other)
		return Set(other);

	return (*this);
}

Mat4f &Mat4f::operator=(Mat3f &matrix3x3)
{
	return Set(matrix3x3);
}


bool Mat4f::operator==(const Mat4f &other) const
{
	int equalities = 0;

	for (int i = 0; i < SIZE_4; i++)
	{
		for (int j = 0; j < SIZE_4; j++)
		{
			equalities += (matrix[i][j] == other.matrix[i][j]) ? 1 : 0;
		}
	}

	return equalities == (SIZE_4 * SIZE_4);
}

bool Mat4f::operator!=(const Mat4f & other) const
{
	return !(*this == other);
}

Mat4f Mat4f::operator+(const Mat4f &other) const
{
	return (Mat4f(*this) += other);
}

Mat4f Mat4f::operator+(float _matrix[SIZE_4][SIZE_4]) const
{
	return (Mat4f(*this) += _matrix);
}

Mat4f Mat4f::operator+(float scalar) const
{
	return (Mat4f(*this) += scalar);
}

Mat4f &Mat4f::operator+=(const Mat4f &other)
{
	matrix[0][0] += other.matrix[0][0];
	matrix[0][1] += other.matrix[0][1];
	matrix[0][2] += other.matrix[0][2];
	matrix[0][3] += other.matrix[0][3];

	matrix[1][0] += other.matrix[1][0];
	matrix[1][1] += other.matrix[1][1];
	matrix[1][2] += other.matrix[1][2];
	matrix[1][3] += other.matrix[1][3];

	matrix[2][0] += other.matrix[2][0];
	matrix[2][1] += other.matrix[2][1];
	matrix[2][2] += other.matrix[2][2];
	matrix[2][3] += other.matrix[2][3];

	matrix[3][0] += other.matrix[3][0];
	matrix[3][1] += other.matrix[3][1];
	matrix[3][2] += other.matrix[3][2];
	matrix[3][3] += other.matrix[3][3];

	return (*this);
}

Mat4f &Mat4f::operator+=(float _matrix[SIZE_4][SIZE_4])
{
	matrix[0][0] += _matrix[0][0];
	matrix[0][1] += _matrix[0][1];
	matrix[0][2] += _matrix[0][2];
	matrix[0][3] += _matrix[0][3];

	matrix[1][0] += _matrix[1][0];
	matrix[1][1] += _matrix[1][1];
	matrix[1][2] += _matrix[1][2];
	matrix[1][3] += _matrix[1][3];

	matrix[2][0] += _matrix[2][0];
	matrix[2][1] += _matrix[2][1];
	matrix[2][2] += _matrix[2][2];
	matrix[2][3] += _matrix[2][3];

	matrix[3][0] += _matrix[3][0];
	matrix[3][1] += _matrix[3][1];
	matrix[3][2] += _matrix[3][2];
	matrix[3][3] += _matrix[3][3];

	return (*this);
}

Mat4f &Mat4f::operator+=(float scalar)
{
	matrix[0][0] += scalar;
	matrix[0][1] += scalar;
	matrix[0][2] += scalar;
	matrix[0][3] += scalar;
					
	matrix[1][0] += scalar;
	matrix[1][1] += scalar;
	matrix[1][2] += scalar;
	matrix[1][3] += scalar;
					
	matrix[2][0] += scalar;
	matrix[2][1] += scalar;
	matrix[2][2] += scalar;
	matrix[2][3] += scalar;
					
	matrix[3][0] += scalar;
	matrix[3][1] += scalar;
	matrix[3][2] += scalar;
	matrix[3][3] += scalar;

	return (*this);

}

Mat4f Mat4f::operator-(const Mat4f &other) const
{
	return (Mat4f(*this) -= other);
}

Mat4f Mat4f::operator-(float matrix[SIZE_4][SIZE_4]) const
{
	return (Mat4f(*this) -= matrix);
}

Mat4f Mat4f::operator-(float scalar) const
{
	return (Mat4f(*this) -= scalar);
}

Mat4f &Mat4f::operator-=(const Mat4f &other)
{
	matrix[0][0] -= other.matrix[0][0];
	matrix[0][1] -= other.matrix[0][1];
	matrix[0][2] -= other.matrix[0][2];
	matrix[0][3] -= other.matrix[0][3];
				 
	matrix[1][0] -= other.matrix[1][0];
	matrix[1][1] -= other.matrix[1][1];
	matrix[1][2] -= other.matrix[1][2];
	matrix[1][3] -= other.matrix[1][3];
				 
	matrix[2][0] -= other.matrix[2][0];
	matrix[2][1] -= other.matrix[2][1];
	matrix[2][2] -= other.matrix[2][2];
	matrix[2][3] -= other.matrix[2][3];
				 
	matrix[3][0] -= other.matrix[3][0];
	matrix[3][1] -= other.matrix[3][1];
	matrix[3][2] -= other.matrix[3][2];
	matrix[3][3] -= other.matrix[3][3];

	return (*this);
}

Mat4f &Mat4f::operator-=(float _matrix[SIZE_4][SIZE_4])
{
	matrix[0][0] -= _matrix[0][0];
	matrix[0][1] -= _matrix[0][1];
	matrix[0][2] -= _matrix[0][2];
	matrix[0][3] -= _matrix[0][3];
				 
	matrix[1][0] -= _matrix[1][0];
	matrix[1][1] -= _matrix[1][1];
	matrix[1][2] -= _matrix[1][2];
	matrix[1][3] -= _matrix[1][3];
				 
	matrix[2][0] -= _matrix[2][0];
	matrix[2][1] -= _matrix[2][1];
	matrix[2][2] -= _matrix[2][2];
	matrix[2][3] -= _matrix[2][3];
				 
	matrix[3][0] -= _matrix[3][0];
	matrix[3][1] -= _matrix[3][1];
	matrix[3][2] -= _matrix[3][2];
	matrix[3][3] -= _matrix[3][3];

	return (*this);
}

Mat4f &Mat4f::operator-=(float scalar)
{
	matrix[0][0] -= scalar;
	matrix[0][1] -= scalar;
	matrix[0][2] -= scalar;
	matrix[0][3] -= scalar;
				 
	matrix[1][0] -= scalar;
	matrix[1][1] -= scalar;
	matrix[1][2] -= scalar;
	matrix[1][3] -= scalar;
				 
	matrix[2][0] -= scalar;
	matrix[2][1] -= scalar;
	matrix[2][2] -= scalar;
	matrix[2][3] -= scalar;
				 
	matrix[3][0] -= scalar;
	matrix[3][1] -= scalar;
	matrix[3][2] -= scalar;
	matrix[3][3] -= scalar;

	return (*this);
}

Mat4f Mat4f::operator*(const Mat4f &other) const
{
	return (Mat4f(*this) *= other);
}

Mat4f Mat4f::operator*(Mat3f &matrix3x3) const
{
	return (Mat4f(*this) *= matrix3x3);
}

Mat4f Mat4f::operator*(float _matrix[SIZE_4][SIZE_4]) const
{
	return (Mat4f(*this) *= _matrix);
}

Vec4f Mat4f::operator*(const Vec4f &vector4) const
{
	Vec4f vet(0, 0, 0, 0);

	for (int i = 0; i < SIZE_4; i++)
	{
		for (int j = 0; j < SIZE_4; j++)
		{
			vet[i] += matrix[i][j] * vector4[j];
		}
	}

	return vet;
}

Mat4f Mat4f::operator*(float scalar) const
{
	return (Mat4f(*this) *= scalar);
}

Mat4f &Mat4f::operator*=(const Mat4f &other)
{
	Mat4f aux(*this);

	for (int count = 0; count < SIZE_4; count++)
	{
		for (int i = 0; i < SIZE_4; i++)
		{
			matrix[count][i] = 0;
			for (int j = 0; j < SIZE_4; j++)
			{
				matrix[count][i] += aux.matrix[count][j] * other.matrix[j][i];
			}
		}
	}

	return (*this);
}

Mat4f &Mat4f::operator*=(Mat3f &matrix)
{
	Mat4f aux(*this);

	for (int count = 0; count < SIZE_3; count++)
	{
		for (int i = 0; i < SIZE_3; i++)
		{
			matrix[count][i] = 0;
			for (int j = 0; j < SIZE_3; j++)
			{
				matrix[count][i] = aux[count][j] * matrix[j][i];
			}
		}
	}

	return (*this);
}

Mat4f &Mat4f::operator*=(float _matrix[SIZE_4][SIZE_4])
{
	Mat4f aux(*this);

	for (int count = 0; count < SIZE_4; count++)
	{
		for (int i = 0; i < SIZE_4; i++)
		{
			matrix[count][i] = 0;
			for (int j = 0; j < SIZE_4; j++)
			{
				matrix[count][i] += aux.matrix[count][j] * _matrix[j][i];
			}
		}
	}

	return (*this);
}

Mat4f &Mat4f::operator*=(float scalar)
{
	for (int i = 0; i < SIZE_4; i++)
	{
		for (int j = 0; j < SIZE_4; j++)
		{
			matrix[i][j] *= scalar;
		}
	}

	return (*this);
}

Mat4f Mat4f::operator/(float scalar) const
{
	return (Mat4f(*this) /= scalar);
}

Mat4f &Mat4f::operator/=(float scalar)
{
	matrix[0][0] /= scalar;
	matrix[0][1] /= scalar;
	matrix[0][2] /= scalar;
	matrix[0][3] /= scalar;
				 
	matrix[1][0] /= scalar;
	matrix[1][1] /= scalar;
	matrix[1][2] /= scalar;
	matrix[1][3] /= scalar;
				 
	matrix[2][0] /= scalar;
	matrix[2][1] /= scalar;
	matrix[2][2] /= scalar;
	matrix[2][3] /= scalar;
				 
	matrix[3][0] /= scalar;
	matrix[3][1] /= scalar;
	matrix[3][2] /= scalar;
	matrix[3][3] /= scalar;

	return (*this);
}

Mat4f &Mat4f::Scale(const Vec3f &vector)
{
	for (int i = 0; i < SIZE_3; i++)
	{
		for (int j = 0; j < SIZE_4; j++)
		{
			matrix[i][j] *= vector[i];
		}
	}

	return (*this);
}

Mat4f &Mat4f::Scale(float _x, float _y, float _z)
{
	return Scale(Vec3f(_x, _y, _z));
}

Mat4f &Mat4f::Scale(float size)
{
	return Scale(Vec3f(size, size, size));
}

Mat4f &Mat4f::RotateX(float angle)
{
	float Cos = cos(angle);
	float Sin = sin(angle);

	Mat4f aux(1,  0,   0,   0,
			   0,  Cos, Sin, 0,
			   0, -Sin, Cos, 0,
			   0,  0,   0,   1);

	aux *= *this;
	return (*this = aux);
}

Mat4f &Mat4f::RotateY(float p_angle)
{
	float Cos = cos(p_angle);
	float Sin = sin(p_angle);

	Mat4f aux(Cos, 0, -Sin, 0,
			   0,   1,  0,   0,
			   Sin, 0,  Cos, 0,
			   0,   0,  0,   1);

	aux *= *this;
	return (*this = aux);
}

Mat4f &Mat4f::RotateZ(float angle)
{
	float Cos = cos(angle);
	float Sin = sin(angle);

	Mat4f aux(Cos, Sin, 0, 0,
			  -Sin, Cos, 0, 0,
			   0,   0,   1, 0,
			   0,   0,   0, 1);

	aux *= *this;
	return (*this = aux);
}

Mat4f &Mat4f::RotateXYZ(float angleX, float angleY, float angleZ)
{
	return RotateX(angleX).RotateY(angleY).RotateZ(angleZ);
}

Mat4f &Mat4f::Translate(const Vec3f &distance)
{
	Vec4f vet(0, 0, 0, 0);

	for (int i = 0; i < SIZE_4; i++)
	{
		for (int j = 0; j < SIZE_4; j++)
		{
			if (j == SIZE_3)
			{
				vet[i] += matrix[j][i];
				continue;
			}
			vet[i] += distance[j] * matrix[j][i];
		}
	}

	return Set(vet);
}

Mat4f &Mat4f::Translate(float _x, float _y, float _z)
{
	return Translate(Vec3f(_x, _y, _z));
}

Mat4f &Mat4f::Ortho(const float left, const float right, const float bottom, 
					  const float top,  const float near,  const float far)
{
	Vec3f aux1(2.0f / (right - left  ),
			    2.0f / (top   - bottom),
			   -2.0f / (far   - near  ) );
	Vec3f aux2(-(right + left  ) / (right - left  ),
				-(top   + bottom) / (top   - bottom),
				-(far   + near  ) / (far   - near  ) );

	for (int i = 0; i < SIZE_4; i++)
	{
		for (int j = 0; j < SIZE_3; j++)
		{
			matrix[3][i] += matrix[j][i] * aux2[j];
		}
	}

	for (int i = 0; i < SIZE_3; i++)
	{
		for (int j = 0; j < SIZE_4; j++)
		{
			matrix[i][j] *= aux1[i];
		}
	}

	return (*this);
}

Mat4f &Mat4f::Perspective(const Vec4f &parameters)
{
	float h = tanf(parameters.x * 0.5f) * parameters.z;
	float w = h *  parameters.y;

	Vec4f aux(parameters.z/w,
			   parameters.z/h,
			  -(parameters.w + parameters.z)/
			   (parameters.w - parameters.z),
			  -2.0f * parameters.w * parameters.z/
			         (parameters.w - parameters.z));

	Mat4f auxM;

	for (int i = 0; i < SIZE_4; i++)
	{
		auxM[2][i] = matrix[2][i] * aux[2] - matrix[3][i];
	}

	for (int i = 0; i < SIZE_4; i++)
	{
		for (int j = 0; j < SIZE_4; j++)
		{
			if (i == 2)
				i++;
			if (i == 3)
			{
				auxM[i][j] = matrix[i-1][j] * aux[i];
				continue;
			}

			auxM[i][j] = matrix[i][j] * aux[i];
		}
	}

	return (*this = auxM);
}

Mat4f &Mat4f::Perspective(const float fovy, const float aspect, 
							const float near, const float far)
{
	return Perspective(Vec4f(fovy, aspect, near, far));
}

Mat4f &Mat4f::LookAt(const Vec3f &eye, const Vec3f &center, const Vec3f &up)
{
	Vec3f dir = center - eye;
	dir.Normalize();

	//float length = 1.0f / (eye - center).Size();// p_center->GetModuleDistance(p_eye);
	//dir *= length;

	Vec3f right(dir.y * up.z - dir.z * up.y,
				dir.z * up.x - dir.x * up.z,
				dir.x * up.y - dir.y * up.x);
	right.Normalize();

	Vec3f resultUp(right.y * dir.z - right.z * dir.y,
				   right.z * dir.x - right.x * dir.z,
				   right.x * dir.y - right.y * dir.x);

	Mat3f auxM(right.x, resultUp.x, -dir.x,
				right.y, resultUp.y, -dir.y,
				right.z, resultUp.z, -dir.z);

	Vec3f auxV(-right.x    * eye.x      - right.y   * 
				 eye.y      - right.z    * eye.z       ,
				-resultUp.x * eye.x      - resultUp.y * 
				 eye.y      - resultUp.z * eye.z       ,
				 dir.x      * eye.x      + dir.y      * 
				 eye.y      + dir.z      * eye.z        );


	Mat4f resulM;
	resulM.Null();

	for (int c = 0; c < SIZE_3; c++)
	{
		for (int i = 0; i < SIZE_4; i++)
		{
			for (int j = 0; j < SIZE_3; j++)
				resulM.matrix[c][i] += auxM[c][j] * matrix[j][i];
		}
	}
	for (int i = 0; i < SIZE_4; i++)
	{
		for (int j = 0; j < SIZE_3; j++)
			resulM.matrix[3][i] += matrix[j][i] * auxV[j];

		resulM.matrix[3][i] += matrix[3][i];
	}

	*this = resulM;

	return (*this);
}

Mat4f &Mat4f::LookAt(const float xEye,    const float yEye,    const float zEye, 
					   const float xCenter, const float yCenter, const float zCenter, 
					   const float xUp,     const float yUp,     const float zUp)
{
	return LookAt(Vec3f(xEye,    yEye,    zEye   ),
				  Vec3f(xCenter, yCenter, zCenter),
				  Vec3f(xUp,     yUp,     zUp    ) );
}

Mat4f RT::operator*(float scalar, const Mat4f &matrix)
{
	return (matrix * scalar);
}

Mat4f RT::operator+(float scalar, const Mat4f &matrix)
{
	return (matrix + scalar);
}

Mat4f RT::operator-(float scalar, const Mat4f &matrix)
{
	return (matrix - scalar);
}
