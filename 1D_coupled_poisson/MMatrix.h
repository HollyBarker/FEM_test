#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H

#include <vector>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<fstream>
#include<string>

// Class that represents a mathematical matrix
class MMatrix
{

friend std::ostream& operator<< (std::ostream& os, MMatrix A)
{
	for (int i=0;i<A.Rows();i++)
	{
		for (int j=0;j<A.Cols();j++)
		{
			os.precision(8);
			os.width(15);os<<A(i,j);
		}
		os<<std::endl;
	}
	return os;
}

public:
	// constructors
	MMatrix() : nRows(0), nCols(0) {}
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x) {}

	// set all matrix entries equal to a double
	MMatrix &operator=(double x)
	{
		for (int i = 0; i < nRows * nCols; i++) A[i] = x;
		return *this;
	}

	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const
	{
		return A[j + i * nCols];
	}

	// access element, indexed by (row, column) [lvalue]
	double &operator()(int i, int j)
	{
		return A[j + i * nCols];
	}

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }

	//set size
	void setRows(int n){nRows=n;}
	void setCols(int n){nCols=n;}


private:
	unsigned int nRows, nCols;
	std::vector<double> A;
};
inline MVector operator*(const MMatrix& A, const MVector& x)
{
	MVector b(A.Rows());
	if (A.Cols() == x.size())
	{
		double sum;
		for (int i=0;i<A.Rows();i++)
		{
			sum=0;
			for (int j=0;j<A.Cols();j++)
			{
				sum+=A(i,j)*x[j];
			}
			b[i]=sum;
		}
		return b;
	}
	else 
	{
		std::cout<<"ERROR: Attempted matrix * vector multiplication where number of matrix columns != vector length."<<std::endl;
		exit(1);
	}
}
#endif

