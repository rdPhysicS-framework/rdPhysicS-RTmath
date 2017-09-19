#include "RTmat3.h"
#include "RTmat4.h"
#include<iostream>
#include"Util.h"

using namespace RT;

Mat3f::Mat3f(float m0x0, float m0x1, float m0x2, 
			   float m1x0, float m1x1, float m1x2, 
			   float m2x0, float m2x1, float m2x2)
{
	Set(m0x0, m0x1, m0x2,
		m1x0, m1x1, m1x2,
		m2x0, m2x1, m2x2);
}

Mat3f::Mat3f(const Mat3f &matrix3x3)
{
	memcpy(matrix, (void*)matrix3x3.matrix, sizeof(matrix));
}

Mat3f::Mat3f(float _matrix[SIZE_3]
						    [SIZE_3]) 
{
	memcpy(matrix, (void*)_matrix, sizeof(matrix));
}

Mat3f::~Mat3f()
{}

Mat3f &Mat3f::Set( float m0x0, float m0x1, float m0x2, 
				   float m1x0, float m1x1, float m1x2, 
				   float m2x0, float m2x1, float m2x2)
{
	matrix[0][0] = m0x0; matrix[0][1] = m0x1; matrix[0][2] = m0x2;
	matrix[1][0] = m1x0; matrix[1][1] = m1x1; matrix[1][2] = m1x2;
	matrix[2][0] = m2x0; matrix[2][1] = m2x1; matrix[2][2] = m2x2;

	return (*this);
}

Mat3f &Mat3f::Set(const Mat3f &other)
{
	memcpy(matrix, (void*)other.matrix, sizeof(matrix));

	return (*this);
}

Mat3f &Mat3f::Set(float matrix3x3[SIZE_3]
								 [SIZE_3])
{
	memcpy(matrix, (void*)matrix3x3, sizeof(matrix));

	return (*this);
}

Mat3f &Mat3f::Set(const Vec3f &vector)
{
	matrix[0][2] = vector.x;
	matrix[1][2] = vector.y;
	matrix[2][2] = vector.z;

	return (*this);
}

Mat3f &RT::Mat3f::Set(const Vec2f & vector)
{
	matrix[0][2] = vector.x;
	matrix[1][2] = vector.y;

	return (*this);
}

float Mat3f::Getixj(int i, int j) const
{
	return matrix[i][j];
}

int Mat3f::GetWidth() const
{
	return SIZE_3;
}

int Mat3f::GetHeight() const
{
	return SIZE_3;
}

int Mat3f::Length() const
{
	return SIZE_3 * SIZE_3;
}

Mat3f &Mat3f::Transpose()
{
	Mat3f aux(*this);

	matrix[0][0] = aux[0][0];
	matrix[0][1] = aux[1][0];
	matrix[0][2] = aux[2][0];
	matrix[1][0] = aux[0][1];
	matrix[1][1] = aux[1][1];
	matrix[1][2] = aux[2][1];
	matrix[2][0] = aux[0][2];
	matrix[2][1] = aux[1][2];
	matrix[2][2] = aux[2][2];

	return (*this);
}

Mat3f &Mat3f::Inverse()
{
	float s = 1 / Determinant();

	return Set((matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) * s,
			   (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) * s,
			   (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * s,
			   														 
			   (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) * s,
			   (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * s,
			   (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) * s,
			   														 
			   (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) * s,
			   (matrix[0][2] * matrix[2][0] - matrix[0][0] * matrix[2][1]) * s,
			   (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * s );
}

Mat3f &Mat3f::Identity()
{
	return Set(1.0f, 0.0f, 0.0f,
			   0.0f, 1.0f, 0.0f,
		       0.0f, 0.0f, 1.0f);
}

Mat3f &Mat3f::Opposed()
{
	return (*this * -1);
}

Mat3f &Mat3f::Null()
{
	matrix[0][0] = matrix[0][1] = matrix[0][2] = 0;
	matrix[1][0] = matrix[1][1] = matrix[1][2] = 0;
	matrix[2][0] = matrix[2][1] = matrix[2][2] = 0;

	return (*this);
}

float Mat3f::Determinant() const
{
	return matrix[0][0] *
		 ((matrix[1][1] * matrix[2][2]) -
		 	(matrix[1][2] * matrix[2][1])) +
		 
		 matrix[1][0] *
		 ((matrix[0][1] * matrix[2][2]) -
		 	(matrix[0][2] * matrix[2][1])) +
		 
		 matrix[2][0] *
		 ((matrix[0][1] * matrix[1][2]) -
		 	(matrix[0][2] * matrix[1][1]));
}

Mat3f &Mat3f::operator=(const Mat3f &other)
{
	return Set(other);
}

bool Mat3f::operator==(const Mat3f &other) const
{
	int equalities = 0;

	for (int i = 0; i < SIZE_3; i++)
	{
		for (int j = 0; j < SIZE_3; j++)
		{
			equalities += (matrix[i][j] == other.matrix[i][j]) ? 1 : 0;
		}
	}

	return equalities == (SIZE_3 * SIZE_3);
}

bool Mat3f::operator!=(const Mat3f & other) const
{
	return !(*this == other);
}

Mat3f &Mat3f::operator+=(const Mat3f &other)
{
	matrix[0][0] += other.matrix[0][0];
	matrix[0][1] += other.matrix[0][1];
	matrix[0][2] += other.matrix[0][2];

	matrix[1][0] += other.matrix[1][0];
	matrix[1][1] += other.matrix[1][1];
	matrix[1][2] += other.matrix[1][2];

	matrix[2][0] += other.matrix[2][0];
	matrix[2][1] += other.matrix[2][1];
	matrix[2][2] += other.matrix[2][2];

	return (*this);
}

Mat3f &Mat3f::operator+=(float _matrix[SIZE_3]
									  [SIZE_3])
{
	matrix[0][0] += _matrix[0][0];
	matrix[0][1] += _matrix[0][1];
	matrix[0][2] += _matrix[0][2];
					
	matrix[1][0] += _matrix[1][0];
	matrix[1][1] += _matrix[1][1];
	matrix[1][2] += _matrix[1][2];
					
	matrix[2][0] += _matrix[2][0];
	matrix[2][1] += _matrix[2][1];
	matrix[2][2] += _matrix[2][2];

	return (*this);
}

Mat3f & RT::Mat3f::operator+=(float scalar)
{
	matrix[0][0] += scalar;
	matrix[0][1] += scalar;
	matrix[0][2] += scalar;
					
	matrix[1][0] += scalar;
	matrix[1][1] += scalar;
	matrix[1][2] += scalar;
					
	matrix[2][0] += scalar;
	matrix[2][1] += scalar;
	matrix[2][2] += scalar;

	return (*this);
}

Mat3f Mat3f::operator+(const Mat3f &other) const
{
	return (Mat3f(*this) += other);
}

Mat3f Mat3f::operator+(float _matrix[SIZE_3]
									  [SIZE_3]) const
{
	return (Mat3f(*this) += _matrix);
}

Mat3f RT::Mat3f::operator+(float scalar) const
{
	return (Mat3f(*this) += scalar);
}

Mat3f &Mat3f::operator-=(const Mat3f &other)
{
	matrix[0][0] -= other.matrix[0][0];
	matrix[0][1] -= other.matrix[0][1];
	matrix[0][2] -= other.matrix[0][2];
				 
	matrix[1][0] -= other.matrix[1][0];
	matrix[1][1] -= other.matrix[1][1];
	matrix[1][2] -= other.matrix[1][2];
				 
	matrix[2][0] -= other.matrix[2][0];
	matrix[2][1] -= other.matrix[2][1];
	matrix[2][2] -= other.matrix[2][2];

	return (*this);
}

Mat3f &Mat3f::operator-=(float _matrix[SIZE_3]
										[SIZE_3])
{
	matrix[0][0] -= _matrix[0][0];
	matrix[0][1] -= _matrix[0][1];
	matrix[0][2] -= _matrix[0][2];
				 
	matrix[1][0] -= _matrix[1][0];
	matrix[1][1] -= _matrix[1][1];
	matrix[1][2] -= _matrix[1][2];
				 
	matrix[2][0] -= _matrix[2][0];
	matrix[2][1] -= _matrix[2][1];
	matrix[2][2] -= _matrix[2][2];

	return (*this);
}

Mat3f &RT::Mat3f::operator-=(float scalar)
{
	matrix[0][0] -= scalar;
	matrix[0][1] -= scalar;
	matrix[0][2] -= scalar;
				 
	matrix[1][0] -= scalar;
	matrix[1][1] -= scalar;
	matrix[1][2] -= scalar;
				 
	matrix[2][0] -= scalar;
	matrix[2][1] -= scalar;
	matrix[2][2] -= scalar;

	return (*this);
}

Mat3f Mat3f::operator-(const Mat3f &other) const
{
	return (Mat3f(*this) -= other);
}

Mat3f Mat3f::operator-(float _matrix[SIZE_3]
									  [SIZE_3]) const
{
	return (Mat3f(*this) -= _matrix);
}

Mat3f RT::Mat3f::operator-(float scalar) const
{
	return (Mat3f(*this) -= scalar);
}

Mat3f &Mat3f::operator*=(const Mat3f &other)
{
	Mat3f aux(*this);

	for (int count = 0; count < SIZE_3; count++)
	{
		for (int i = 0; i < SIZE_3; i++)
		{
			matrix[count][i] = 0;
			for (int j = 0; j < SIZE_3; j++)
			{
				matrix[count][i] += aux[count][j] * other.matrix[j][i];
			}
		}
	}

	return (*this);
}

Mat3f &Mat3f::operator*=(Mat4f &other)
{
	Mat3f aux(*this);

	for (int count = 0; count < SIZE_3; count++)
	{
		for (int i = 0; i < SIZE_3; i++)
		{
			matrix[count][i] = 0;
			for (int j = 0; j < SIZE_3; j++)
			{
				matrix[count][i] += aux[count][j] * other[j][i];
			}
		}
	}

	return (*this);
}

Mat3f &Mat3f::operator*=(float _matrix[SIZE_3]
										[SIZE_3])
{
	return *this *= _matrix;
}

Mat3f Mat3f::operator*(const Mat3f &other) const
{
	return (Mat3f(*this) *= other);
}

Mat3f Mat3f::operator*(Mat4f &other) const
{
	return (Mat3f(*this) *= other);
}

Mat3f Mat3f::operator*(float _matrix[SIZE_3]
									[SIZE_3]) const
{
	return (Mat3f(*this) *= _matrix);
}

Vec3f Mat3f::operator*(const Vec3f &vector) const
{
	Vec3f vec(0.0f, 0.0f, 0.0f);

	for (int i = 0; i < SIZE_3; i++)
	{
		for (int j = 0; j < SIZE_3; j++)
		{
			vec[i] += (matrix[i][j] * vector[j]);
		}
	}

	return vec;
}

Mat3f &Mat3f::operator*=(float scalar)
{
	matrix[0][0] *= scalar;
	matrix[0][1] *= scalar;
	matrix[0][2] *= scalar;
				 	
	matrix[1][0] *= scalar;
	matrix[1][1] *= scalar;
	matrix[1][2] *= scalar;
				 	
	matrix[2][0] *= scalar;
	matrix[2][1] *= scalar;
	matrix[2][2] *= scalar;

	return (*this);
}

Mat3f Mat3f::operator*(float scalar) const
{
	return (Mat3f(*this) *= scalar);
}

Mat3f &Mat3f::operator/=(float scalar)
{
	matrix[0][0] /= scalar;
	matrix[0][1] /= scalar;
	matrix[0][2] /= scalar;
				
	matrix[1][0] /= scalar;
	matrix[1][1] /= scalar;
	matrix[1][2] /= scalar;
				
	matrix[2][0] /= scalar;
	matrix[2][1] /= scalar;
	matrix[2][2] /= scalar;

	return (*this);
}

Mat3f Mat3f::operator/(float scalar) const
{
	return (Mat3f(*this) /= scalar);
}

Mat3f &Mat3f::Translate(float _x, float _y)
{
	Vec3f aux = *this * Vec3f(_x, _y, 1.f);
	return Set(aux);
}

Mat3f &Mat3f::Translate(const Vec2f &vector)
{
	return Translate(vector.x, vector.y);
}

Mat3f &Mat3f::Scale(float scale)
{
	return Scale(scale, scale);
}

Mat3f &Mat3f::Scale(float scalex, float scaley)
{
	return Scale(Vec2f(scalex, scaley));
}

Mat3f &Mat3f::Scale(const Vec2f &scales)
{
	for (int i = 0; i < SIZE_2; i++)
	{
		for (int j = 0; j < SIZE_2; j++)
		{
			matrix[i][j] *= scales[i];
		}
	}

	return (*this);
}

Mat3f &Mat3f::Rotate(float angle)
{
	float Cos = cos(angle);
	float Sin = sin(angle);

	Mat3f aux(Cos, -Sin, 0,
			  Sin,  Cos, 0,
			  0,    0,   1);

	aux *= *this;

	return (*this = aux);
}

Mat3f RT::operator*(float scalar, const Mat3f &matrix)
{
	return (matrix * scalar);
}

Mat3f RT::operator+(float scalar, const Mat3f &matrix)
{
	return (matrix + scalar);
}

Mat3f RT::operator-(float scalar, const Mat3f &matrix)
{
	return (matrix - scalar);
}

std::ostream &RT::operator<<(std::ostream &output, Mat3f &mat)
{
	for (int i = 0; i < SIZE_3; i++)
	{
		for (int j = 0; j < SIZE_3; j++)
			output << mat[i][j] << " ";
		output << std::endl;
	}

	std::cout << std::endl;
	return output;
}

