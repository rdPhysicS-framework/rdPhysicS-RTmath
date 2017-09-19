#include "RTmatrix.h"
#include<iostream>

RTmatrix::RTmatrix(const unsigned int p_width,
			   const unsigned int p_height) :
	    matrix(nullptr), m_width(p_width), m_height(p_height)
{
	Create();
}

RTmatrix::RTmatrix(RTmatrix *p_matrix) : matrix(nullptr)
{
	Create(p_matrix);
}

RTmatrix::RTmatrix(float **p_matrix, 
			   const unsigned int p_width, 
			   const unsigned int p_height) :
		matrix(p_matrix), m_width(p_width), m_height(p_height)
{}

RTmatrix::~RTmatrix()
{
	Destroy();
}

void RTmatrix::Destroy()
{
	for (int i = 0; i < m_width; i++)
	{
		delete matrix[i];
		matrix[i] = nullptr;
	}

	delete[]matrix;
	matrix = nullptr;
}

int RTmatrix::GetWidth() const
{
	return m_width;
}

int RTmatrix::GetHeight() const
{
	return m_height;
}

RTmatrix *RTmatrix::Set(RTmatrix *p_other)
{
	if (p_other)
	{
		for (int i = 0; i < m_height; i++)
		{
			for (int j = 0; j < m_width; j++)
				matrix[i][j] = *((*p_other + i) + j);
		}
	}

	return this;
}

RTmatrix *RTmatrix::Set(float **p_matrix)
{
	if (p_matrix)
	{
		for (int i = 0; i < m_height; i++)
		{
			for (int j = 0; j < m_width; j++)
				matrix[i][j] = p_matrix[i][j];
		}
	}

	return this;
}

RTmatrix *RTmatrix::Set(float *p_matrix)
{
	if (p_matrix)
	{
		int count = 0;

		for (int i = 0; i < m_height; i++)
		{
			for (int j = 0; j < m_width; j++)
			{
				matrix[i][j] = p_matrix[count];
				count++;
			}
		}
	}

	return this;
}

float **RTmatrix::GetMatrix() const
{
	return matrix;
}

RTmatrix *RTmatrix::Create()
{
	if (!matrix)
	{
		matrix = new float*[m_height];

		for (int i = 0; i < m_height; i++)
		{
			matrix[i] = new float[m_width];
		}
	}

	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			if (j == 0)
			{
				matrix[i][j] = 1;
				continue;
			}

			matrix[i][j] = 0;
		}
	}

	return this;
}

RTmatrix *RTmatrix::Create(RTmatrix *p_matrix)
{
	if(!p_matrix->GetMatrix())
		return this;

	m_width = p_matrix->GetWidth();
	m_height = p_matrix->GetHeight();

	if (!matrix)
	{
		matrix = new float*[m_height];

		for (int i = 0; i < m_height; i++)
		{
			matrix[i] = new float[m_width];
		}
	} 


	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
			matrix[i][j] = *((*p_matrix + i) + j);
	}

	return this;
}

RTmatrix *RTmatrix::Create(float **p_matrix)
{
	return this;
}

float *RTmatrix::operator+(int p_index)
{
	return matrix[p_index];
}

float *RTmatrix::operator[](const unsigned int p_index)
{
	return matrix[p_index];
}

RTmatrix *RTmatrix::operator=(RTmatrix *p_other)
{
	if (this != p_other)
	{
		Set(p_other);
	}

	return this;
}

RTmatrix RTmatrix::operator+(float p_scalar)
{
	return (RTmatrix(this) += p_scalar);
}

RTmatrix RTmatrix::operator+(RTmatrix *p_matrix)
{
	return RTmatrix(this) += p_matrix;
}

RTmatrix *RTmatrix::operator+=(float p_scalar)
{
	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			matrix[i][j] += p_scalar;
		}
	}

	return this;
}

RTmatrix *RTmatrix::operator+=(RTmatrix *p_matrix)
{
	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			matrix[i][j] += *((*p_matrix + i) + j);
		}
	}

	return this;
}

RTmatrix RTmatrix::operator-(float p_scalar)
{
	return RTmatrix(this) -= p_scalar;
}

RTmatrix RTmatrix::operator-(RTmatrix *p_matrix)
{
	return RTmatrix(this) -= p_matrix;
}

RTmatrix *RTmatrix::operator-=(float p_scalar)
{
	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			matrix[i][j] -= p_scalar;
		}
	}

	return this;
}

RTmatrix *RTmatrix::operator-=(RTmatrix *p_matrix)
{
	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			matrix[i][j] -= *((*p_matrix + i) + j);// p_matrix->GetMatrix()[i][j];
		}
	}

	return this;
}

RTmatrix RTmatrix::operator*(float p_scalar)
{
	return RTmatrix(this) *= p_scalar;
}

RTmatrix RTmatrix::operator*(RTmatrix *p_matrix)
{
	if (p_matrix->GetHeight() != m_height ||
		p_matrix->GetWidth() != m_width)
	{
		return this;
	}

	return RTmatrix(this) *= p_matrix;
}

RTmatrix *RTmatrix::operator*=(float p_scalar)
{
	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			matrix[i][j] *= p_scalar;
		}
	}

	return this;
}

RTmatrix *RTmatrix::operator*=(RTmatrix *p_matrix)
{
	if (p_matrix->GetHeight() != m_height ||
		p_matrix->GetWidth() != m_width)
	{
		return this;
	}

	RTmatrix aux(this);


	for (int count = 0; count < m_height; count++)
	{
		for (int i = 0; i < m_width; i++)
		{
			matrix[count][i] = 0;

			for (int j = 0; j < m_width; j++)
			{
				matrix[count][i] += aux[count][j] * (*((*p_matrix + j) + i));
			}
		}
	}

	return this;
}

RTmatrix *RTmatrix::operator/=(float p_scalar)
{
	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
			matrix[i][j] /= p_scalar;
	}

	return this;
}

RTmatrix RTmatrix::operator/(float p_scalar)
{
	return RTmatrix(this) /= p_scalar;
}

RTmatrix *RTmatrix::SetTranspose(RTmatrix *p_dest)
{
	if (this == p_dest)
		return SetTranspose();

	return nullptr;
}

RTmatrix *RTmatrix::SetOpposite(RTmatrix *p_dest)
{
	return nullptr;
}

RTmatrix *RTmatrix::SetInverse(RTmatrix *p_dest)
{
	return nullptr;
}

RTmatrix *RTmatrix::SetTranspose()
{
	RTmatrix aux(this);

	Destroy();

	Create(&aux);

	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			matrix[j][i] = aux[i][j];
		}
	}

	return this;
}

RTmatrix *RTmatrix::SetOpposed()
{
	*this *= -1;

	return this;
}

RTmatrix *RTmatrix::SetInverse()
{
	if (m_height != m_width)
		throw std::exception("ERROR...\n");

	float det = CalcDet();

	if (det == 0)
	{
		std::cout << "..." << std::endl;
		return this;
	}

	Adjoint();

	float div = 1 / det;

	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			matrix[i][j] *= div;
		}

	}

	return this;
}

float RTmatrix::CalcDet()
{
	if (m_height == m_width)
	{
		return GetDeterminant(matrix, m_height);
	}

	throw std::exception("ERROR...\n");
}


float RTmatrix::Determinant2X2(float **p_matrix) const
{
	return (p_matrix[0][0] * p_matrix[1][1]) -
		   (p_matrix[0][1] * p_matrix[1][0]);
}

float **RTmatrix::Cofator(float **p_matrix, int p_size, int p_i, int p_j)
{
	if (p_i >= p_size || p_j >= p_size)
	{
		throw std::exception("ERROR...\n");
	}
	if (p_size > 0)
	{
		float **aux = new float*[p_size - 1];

		for (int i = 0; i < p_size - 1; i++)
		{
			aux[i] = new float[p_size - 1];
		}

		int countI = 0, countJ = 0;

		for (int i = 0; i < p_size; i++)
		{
			if (i != p_i)
			{
				for (int j = 0; j < p_size; j++)
				{
					if (j == p_j)
						continue;

					aux[countI][countJ] = p_matrix[i][j];
					countJ++;
				}
				countI++;
				countJ = 0;
			}
		}

		return aux;
	}

	return p_matrix;
}

RTmatrix *RTmatrix::Adjoint(RTmatrix *p_dest)
{
	RTmatrix aux(this);

	float **newMatrix = nullptr;

	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			newMatrix = Cofator(aux.GetMatrix(), m_height, i, j);

			p_dest->GetMatrix()[i][j] = pow(-1, ((i + 1) + (j + 1))) * GetDeterminant(newMatrix, m_height - 1);
		}
	}

	SetTranspose();

	return p_dest;
}

RTmatrix *RTmatrix::Adjoint()
{
	RTmatrix aux(this);

	float **newMatrix = nullptr;

	for (int i = 0; i < m_height; i++)
	{
		for (int j = 0; j < m_width; j++)
		{
			newMatrix = Cofator(aux.GetMatrix(), m_height, i, j);

			matrix[i][j] = pow(-1, ((i + 1) + (j + 1))) * GetDeterminant(newMatrix, m_height - 1);
		}
	}

	SetTranspose();

	return this;
}

float RTmatrix::GetDeterminant(float **p_matrix, int p_size)
{
	if (p_size == 1)
		return **p_matrix;

	float **newMatrix = nullptr;

	float cof;
	float det = 0;

	for (int i = 0; i < p_size; i++)
	{
		cof = p_matrix[i][0];
		newMatrix = Cofator(p_matrix, p_size, i, 0);
		det += (pow(-1, (i + 1 + 1)) * cof) * GetDeterminant(newMatrix, p_size - 1);
	}

	return det;
}
