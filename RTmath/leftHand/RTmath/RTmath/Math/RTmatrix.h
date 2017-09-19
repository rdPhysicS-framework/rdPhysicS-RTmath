#ifndef __RT_MATRIX_H__
#define __RT_MATRIX_H__

class RTmatrix
{
private:
	float **matrix;
	unsigned int m_width;
	unsigned int m_height;

	RTmatrix *Create();
	RTmatrix *Create(RTmatrix *);
	RTmatrix *Create(float **);

	float Determinant2X2(float **) const;
	float **Cofator(float **, int, int, int);
	float GetDeterminant(float **, int);

public:
	RTmatrix(const unsigned int p_width = 2, const unsigned int p_height = 2);
	RTmatrix(RTmatrix *);
	RTmatrix(float **, const unsigned int, const unsigned int);
	~RTmatrix();

	void Destroy();

	int GetWidth() const;
	int GetHeight() const;

	RTmatrix *Set(RTmatrix *p_other);
	RTmatrix *Set(float **p_matrix);
	RTmatrix *Set(float *p_matrix);

	float **GetMatrix() const;

	float *operator+(int p_index);
	float *operator[](const unsigned int);

	RTmatrix *operator=(RTmatrix *p_other);
	RTmatrix operator+(float);
	RTmatrix operator+(RTmatrix *);
	RTmatrix *operator+=(float);
	RTmatrix *operator+=(RTmatrix *);
	RTmatrix operator-(float);
	RTmatrix operator-(RTmatrix *);
	RTmatrix *operator-=(float);
	RTmatrix *operator-=(RTmatrix *);
	RTmatrix operator*(float);
	RTmatrix operator*(RTmatrix *);
	RTmatrix *operator*=(float);
	RTmatrix *operator*=(RTmatrix *);
	RTmatrix *operator/=(float);
	RTmatrix operator/(float);

	RTmatrix *SetTranspose(RTmatrix *);
	RTmatrix *SetOpposite(RTmatrix *);
	RTmatrix *SetInverse(RTmatrix *);
	RTmatrix *Adjoint(RTmatrix *);

	RTmatrix *SetTranspose();
	RTmatrix *SetOpposed();
	RTmatrix *SetInverse();
	RTmatrix *Adjoint();

	float CalcDet();
};
#endif //__RT_MATRIX_H__

