//----------------------------------------
//  Created by Austin Ladshaw on 1/7/14
//  Copyright (c) 2014
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

#ifndef MACAW_HPP_
#define MACAW_HPP_

#include <stdio.h>				//Line to allow cout functionality
#include <math.h>               //Line added to allow usage of the pow (e, x) function
#include <iostream>				//Line to allow for read/write to the console using cpp functions
#include <fstream>				//Line to allow for read/write to and from .txt files
#include <stdlib.h>				//Line need to convert strings to doubles
#include <vector>				//Line needed to use dynamic arrays called vectors
#include <time.h>				//Line needed to display program runtime
#include <float.h>				//Line to allow use of machine precision constants
#include <string>				//Line to allow use of strings as a data type
#include <exception>            //Line to allow use of try-catch statements
#include "error.h"

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#endif

//Matrix Class
template <class T>
class Matrix
{
public:
	
	//Generalized Matrix Operations
	Matrix(int rows, int columns);
	T& operator()(int i, int j);
	T operator()(int i, int j) const;
	Matrix(const Matrix& M);
	Matrix& operator=(const Matrix& M);
	Matrix();
    ~Matrix();
	void set_size(int i, int j);
    void zeros();
    void edit(int i, int j, T value);
    int rows();
    int columns();
	T determinate();
    T norm();
	T sum();
    T inner_product(const Matrix& x);
	Matrix& cofactor(const Matrix& M);
	Matrix operator+(const Matrix& M);
	Matrix operator-(const Matrix& M);
	Matrix operator*(const T);
	Matrix operator/(const T);
	Matrix operator*(const Matrix& M);
	Matrix& transpose(const Matrix& M);
    Matrix& transpose_multiply(const Matrix& MT, const Matrix& v);
	Matrix& adjoint(const Matrix &M);
	Matrix& inverse(const Matrix &M);
	void Display(const std::string Name);
	
	//Specialized Matrix Operations for 1-D FDM
	Matrix& tridiagonalSolve(const Matrix& A, const Matrix& b);
	Matrix& ladshawSolve(const Matrix& A, const Matrix& d);
	Matrix& tridiagonalFill(const T A, const T B, const T C, bool Spherical);
    Matrix& naturalLaplacian3D(int m);
	Matrix& sphericalBCFill(int node, const T coeff, T variable);
	Matrix& ConstantICFill(const T IC);
	
	//Specialized Matrix Operations for SKIMMER
	Matrix& SolnTransform(const Matrix& A, bool Forward);
	T sphericalAvg(double radius, double dr, double bound, bool Dirichlet);
	T IntegralAvg(double radius, double dr, double bound, bool Dirichlet);
	T IntegralTotal(double dr, double bound, bool Dirichlet);
	
	//Specialized Matrix Operations for 1-D Conservation Laws
	Matrix& tridiagonalVectorFill(const std::vector<T> &A, const std::vector<T> &B, const std::vector<T> &C);
	Matrix& columnVectorFill(const std::vector<T> &A);
	Matrix& columnProjection(const Matrix& b, const Matrix& b_old, const double dt, const double dt_old);
	Matrix& dirichletBCFill(int node, const T coeff, T variable);
    
    //Matrix operations for functions focused on Krylov, GMRES, and PJFNK methods
  	Matrix& diagonalSolve(const Matrix& D, const Matrix& v);
    Matrix& upperTriangularSolve(const Matrix& U, const Matrix& v);
    Matrix& lowerTriangularSolve(const Matrix& L, const Matrix& v);
    Matrix& upperHessenberg2Triangular(Matrix& b);
    Matrix& lowerHessenberg2Triangular(Matrix& b);
    Matrix& upperHessenbergSolve(const Matrix& H, const Matrix& v);
    Matrix& lowerHessenbergSolve(const Matrix& H, const Matrix& v);
    Matrix& columnExtract(int j, const Matrix& M);
    Matrix& rowExtract(int i, const Matrix& M);
    Matrix& columnReplace(int j, const Matrix& v);
    Matrix& rowReplace(int i, const Matrix& v);
    void rowShrink();
    void columnShrink();
    void rowExtend(const Matrix& v);
    void columnExtend(const Matrix& v);
    
	
protected:
    int num_rows;
	int num_cols;
	std::vector<T> Data;
};

//Defined methods for the template class
template <class T>
Matrix<T>::Matrix(int rows, int columns)
:
num_rows(rows),
num_cols(columns),
Data(rows * columns)
{
	
}

//For setting values
template <class T>
T& Matrix<T>::operator()(int i, int j)
{
	if (i>=num_rows || j>=num_cols || i<0 || j<0)
	{
		mError(out_of_bounds);
		return Data[0];
	}
	else
		return Data[(i*num_cols)+j];
}

//For accessing values
template <class T>
T Matrix<T>::operator()(int i, int j) const
{
	if (i>=num_rows || j>=num_cols || i<0 || j<0)
	{
		mError(out_of_bounds);
		return Data[0];
	}
	else
	{
		return Data[(i*num_cols)+j];
	}
}

//For copying a matrix at initialization
template <class T>
Matrix<T>::Matrix(const Matrix& M)
:
num_rows(M.num_rows),
num_cols(M.num_cols),
Data(M.num_rows * M.num_cols)
{
	for (int i=0; i<num_rows; i++)
	{
		for (int j=0; j<num_cols; j++)
		{
      		Data[(i*num_cols)+j] = M(i,j);
		}
	}
}

//For setting one matrix equal to anther of same size
template <class T>
Matrix<T> &Matrix<T>::operator=(const Matrix &M)
{
	if (this == &M)
	{
		return *this;
	}
	else
	{
		if (num_rows!=M.num_rows || num_cols!=M.num_cols)
		{
			Data.clear();
			num_rows=M.num_rows;
			num_cols=M.num_cols;
			Data.resize(num_rows*num_cols);
		}
		for (int i=0; i<num_rows; i++)
		{
			for (int j=0; j<num_cols; j++)
			{
				Data[(i*num_cols)+j] = M(i,j);
			}
		}
		return *this;
	}
}

//Default Constructor
template <class T>
Matrix<T>::Matrix()
:
Data(1)
{
	num_cols = 1;
	num_rows = 1;
}

//Default Destructor
template <class T>
Matrix<T>::~Matrix()
{
    Data.clear();
}

//Allocate size
template <class T>
void Matrix<T>::set_size(int i, int j)
{
    if (i <= 0 || j <= 0)
    {
        mError(invalid_size);
        return;
    }
	num_rows = i;
	num_cols = j;
	Data.clear();
	Data.resize(i*j);
}

//Reset all existing matrix entries to zero
template <class T>
void Matrix<T>::zeros()
{
  	for (int n=0; n<this->Data.size(); n++)
		this->Data[n] = 0;
}

//Add/Change an element to a sparse matrix
template <class T>
void Matrix<T>::edit(int i, int j, T value)
{
    if (i>=num_rows || j>=num_cols || i<0 || j<0)
    {
        mError(out_of_bounds);
        return;
    }
    this->Matrix::operator()(i, j) = value;
}

//Function to return the number of rows in the matrix
template <class T>
int Matrix<T>::rows()
{
    return this->num_rows;
}

//Function to return the number of columns in the matrix
template <class T>
int Matrix<T>::columns()
{
    return this->num_cols;
}

//Function to determine the determinate of a matrix
template <class T>
T Matrix<T>::determinate()
{
	T det=0;
	if (num_rows!=num_cols)
	{
		mError(non_square_matrix);
		return 0.0;
	}
	else if (num_rows==2)
	{
		return (Data[0]*Data[3]) - (Data[1]*Data[2]);
	}
	else if (num_rows==1)
	{
		return Data[0];
	}
	else
	{
		int I=0;
		int J=0;
		for(int k=0; k<num_cols; k++)
		{
			//Pivoting
			if (Data[k]==0)
			{
				det = det+0;
				break;
			}
			Matrix<T> temp((num_rows-1), (num_cols-1));
			int r=0;
			int c=0;
			for(int i=0; i<num_rows; i++)
			{
				for(int j=0; j<num_cols; j++)
				{
					if (i!=I && j!=J)
					{
						temp(r,c) = Data[(i*num_cols)+j];
						if (c<(temp.num_cols)-1)
							c++;
					}
				}
				c=0;
				if(i!=I)
				{
					if (r<(temp.num_rows)-1)
						r++;
				}
			}//Filled out temp
			
			det = det + (pow(-1,k)*Data[(0+k)]*temp.determinate());
			if (J<(num_cols-1))
			{
				J++;
			}
			else
			{
				J=0;
				I++;
			}
			temp.Data.clear();
		}//End of column loop
		
		return det;
	}
}

//Calculates the 2-Norm of a matrix or vector
template <class T>
T Matrix<T>::norm()
{
	T norm = 0;
	/*Note: if M is a vector, this returns the 2-Norm,
	 else if M is a matrix it returns the
	 Frobenius Norm of that matrix
	 */
	
	for(int i=0; i<num_rows; i++)
	{
		for (int j=0; j<num_cols; j++)
		{
			norm = norm + pow( fabs(Data[(i*num_cols)+j]) ,2.0);
		}
	}
	norm = sqrt(norm);
	return norm;
}

//Calcuates the sum of all entries
template <class T>
T Matrix<T>::sum()
{
	T sum = 0;
	for(int i=0; i<num_rows; i++)
	{
		for (int j=0; j<num_cols; j++)
		{
			sum = sum + Data[(i*num_cols)+j];
		}
	}
	return sum;
}

//Calculates the inner product of two vectors
template <class T>
T Matrix<T>::inner_product(const Matrix& x)
{
	T prod = 0;
	if (x.num_cols != 1 || num_cols != 1)
	{
		mError(dim_mis_match);
		std::cout << "Inner Product can only be performed on vectors!" << std::endl;
		return prod;
	}
	else if (x.num_rows != num_rows)
	{
		mError(dim_mis_match);
		return prod;
	}
	else
	{
		for (int i=0; i<x.num_rows; i++)
		{
			prod = prod + x(i,0)*Data[(i*num_cols)+0];
		}
	}
	return prod;
}

//Forming the cofactor of a square matrix
template <class T>
Matrix<T> &Matrix<T>::cofactor(const Matrix &M)
{
	if (this == &M)
	{
		mError(arg_matrix_same);
		return *this;
	}
	num_rows = M.num_rows;
	num_cols = M.num_cols;
	this->Data.resize(num_rows*num_cols);
	if (M.num_rows!=M.num_cols)
	{
		mError(non_square_matrix);
		return *this;
	}
	else if (M.num_rows==2)
	{
		this->Data[0] = M.Data[3];
		this->Data[1] = -M.Data[2];
		this->Data[2] = -M.Data[1];
		this->Data[3] = M.Data[0];
		return *this;
	}
	else if (M.num_rows==1)
	{
		return *this;
	}
	else
	{
		int I=0;	//row index of cofactor
		int J=0;	//col index of cofactor
		//Loop for all elements in cofactor
		for (int k=0; k<(M.num_rows * M.num_cols); k++)
		{
			Matrix<T> temp((M.num_rows-1), (M.num_cols-1));
			
			int r=0;	//row index of temp
			int c=0;	//col index of temp
			
			//Loop for all rows of M
			for(int i=0; i<M.num_rows; i++)
			{
				//Loop for all cols of M
				for(int j=0; j<M.num_cols; j++)
				{
					if (i!=I && j!=J)
					{
						temp(r,c) = M.Data[(i*M.num_cols)+j];
						if (c<(temp.num_cols)-1)
							c++;
					}
				}
				c=0;
				if(i!=I)
				{
					if (r<(temp.num_rows)-1)
						r++;
				}
			}
			this->Data[(I*num_cols)+J] = pow(-1,(I+J))*temp.determinate();
			if (J<(M.num_cols-1))
			{
				J++;
			}
			else
			{
				J=0;
				I++;
			}
			temp.Data.clear();
		}
		return *this;
	}
}

//Matrix addition
template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix &M)
{
	if(num_rows!=M.num_rows && num_cols!=M.num_cols)
	{
		mError(dim_mis_match);
		return *this;
	}
	else
	{
		Matrix<T> temp(num_rows,num_cols);
		for (int i=0; i<num_rows; i++)
		{
			for (int j=0; j<num_cols; j++)
			{
				temp.Data[(i*num_cols)+j] = this->Data[(i*num_cols)+j] + M.Data[(i*M.num_cols)+j];
			}
		}
		return temp;
	}
}

//Matrix subtraction
template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix &M)
{
	if(num_rows!=M.num_rows && num_cols!=M.num_cols)
	{
		mError(dim_mis_match);
		return *this;
	}
	else
	{
		Matrix<T> temp(num_rows,num_cols);
		for (int i=0; i<num_rows; i++)
		{
			for (int j=0; j<num_cols; j++)
			{
				temp.Data[(i*num_cols)+j] = this->Data[(i*num_cols)+j] - M.Data[(i*M.num_cols)+j];
			}
		}
		return temp;
	}
}

//Matrix scalar multiplication
template <class T>
Matrix<T> Matrix<T>::operator*(const T a)
{
	Matrix<T> temp(num_rows,num_cols);
	for (int i=0; i<num_rows; i++)
	{
		for (int j=0; j<num_cols; j++)
		{
			temp.Data[(i*num_cols)+j] = a * this->Data[(i*num_cols)+j];
		}
	}
	return temp;
}

//Matrix scalar division
template <class T>
Matrix<T> Matrix<T>::operator/(const T a)
{
	Matrix<T> temp(num_rows,num_cols);
	for (int i=0; i<num_rows; i++)
	{
		for (int j=0; j<num_cols; j++)
		{
			temp.Data[(i*num_cols)+j] = this->Data[(i*num_cols)+j] / a;
		}
	}
	return temp;
}

//Matrix  multiplication
template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix &M)
{
	if (num_cols!=M.num_rows)
	{
		mError(dim_mis_match);
		return *this;
	}
	else
	{
		Matrix<T> temp(num_rows, M.num_cols);
		temp.num_rows = num_rows;
		temp.num_cols = M.num_cols;
		temp.Data.resize(num_rows*M.num_cols);
		
		int j=0;
		for(int i=0; i<num_rows; i++)
		{
			for (int J=0; J<M.num_cols; J++)
			{
				temp.edit(i, J, 0);
				j=0;
				for (int I=0; I<M.num_rows; I++)
				{
					temp.edit(i, J,temp(i,J) + (this->Data[(i*num_cols)+j] * M.Data[(I*M.num_cols)+J]) );
					j++;
				}
			}
		}
		return temp;
	}
}

//Transpose of a matrix
template <class T>
Matrix<T> &Matrix<T>::transpose(const Matrix &M)
{
	if (this == &M)
	{
		mError(arg_matrix_same);
		return *this;
	}
	//A will be indexed by rows and cols
	this->set_size(M.num_cols, M.num_rows);
	
	for (int i=0; i<num_rows; i++)
	{
		for (int j=0; j<num_cols; j++)
		{
			this->Data[(i*num_cols)+j] = M.Data[(j*M.num_cols)+i];
		}
	}
	
	return *this;
}

//Transposes the matrix MT then multiplies by matrix v in one step
template <class T>
Matrix<T> &Matrix<T>::transpose_multiply(const Matrix &MT, const Matrix &v)
{
	if (MT.num_rows != v.num_rows)
	{
		mError(dim_mis_match);
		return *this;
	}
	else if (this == &MT || this == &v)
	{
		mError(arg_matrix_same);
		return *this;
	}
	else
	{
		this->set_size(MT.num_cols, v.num_cols);
		int j=0;
		for(int i=0; i<num_rows; i++)
		{
			for (int J=0; J<v.num_cols; J++)
			{
				this->Data[(i*num_cols)+J] = 0;
				j=0;
				for (int I=0; I<v.num_rows; I++)
				{
					this->Data[(i*num_cols)+J] = this->Data[(i*num_cols)+J] + (MT.Data[(j*MT.num_cols)+i] * v.Data[(I*v.num_cols)+J]);
					j++;
				}
			}
		}
	}
	return *this;
}

//Adjoint of a square matrix
template <class T>
Matrix<T> &Matrix<T>::adjoint(const Matrix &M)
{
	if(num_rows!=num_cols)
	{
		mError(non_square_matrix);
		return *this;
	}
	else if (this == &M)
	{
		mError(arg_matrix_same);
		return *this;
	}
	else
	{
		this->set_size(M.num_rows,M.num_cols);
		this->transpose(this->cofactor(M));
		return *this;
	}
}

//Inverse of a square matrix
template <class T>
Matrix<T> &Matrix<T>::inverse(const Matrix &M)
{
	if(num_rows!=num_cols)
	{
		mError(non_square_matrix);
		return *this;
	}
	else if (this == &M)
	{
		mError(arg_matrix_same);
		return *this;
	}
	else
	{
		this->set_size(M.num_rows,M.num_cols);
		double det_inv;
		Matrix<T> A(M);
		det_inv = 1/A.determinate();
		*this = this->adjoint(A) * det_inv;
		A.Data.clear();
		return *this;
	}
}

//Display Matrix to console
template <class T>
void Matrix<T>::Display(const std::string Name)
{
	std::cout << Name << " = " << std::endl;
	if (num_rows == 0 || num_cols == 0)
	{
		mError(empty_matrix);
		return;
	}
	else
	{
		for (int i=0; i<num_rows; i++)
		{
			std::cout << "\t";
			for (int j=0; j<num_cols; j++)
			{
				std::cout << this->Data[(i*num_cols)+j] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

//Function to solve the tridiagonal matrix using the Thomas Algorithm (only for symmetric matrix)
template <class T>
Matrix<T> &Matrix<T>::tridiagonalSolve(const Matrix& A, const Matrix& b)
{
	if (this == &A || this == &b)
	{
		mError(arg_matrix_same);
		return *this;
	}
	//Solves the system Ax=b when A is a tridiagonal matrix
	this->set_size(b.num_rows, b.num_cols);
	
	//Check dimensions of A and b to ensure they will work out
	if (num_rows!=A.num_cols || num_cols>1)
	{
		mError(dim_mis_match);
		return *this;
	}
	//Check the values in A to ensure it is tridiagonal (Skipped)
	
	/*
	 * solve for x (i.e. this)
	 * a[i] -> 1 to n-1 (lower diag)
	 * c[i] -> 0 to n-2 (upper diag)
	 * b[i] -> 0 to n-1 (diag)
	 * d[i] -> 0 to n-1 (column vector)
	 * cp is c prime and dp is b prime (Ax=b)
	 *
	 * */
	std::vector<T> cp; cp.resize(num_rows-1);
	std::vector<T> dp; dp.resize(num_rows);
	
	//Forward Sweep
	for (int i=0; i<num_rows; i++)
	{
		if (i==0)
		{
			cp[i] = A(i,(i+1))/A(i,i);
			dp[i] = b(i,0)/A(i,i);
		}
		else if (i<(num_rows-1))
		{
			cp[i] = A(i,(i+1)) / (A(i,i) - (cp[(i-1)]*A((i+1),i)) );
			dp[i] = ( b(i,0) - (dp[(i-1)]*A((i+1),i) ) ) / (A(i,i) - (cp[(i-1)]*A((i+1),i)) );
		}
		else
		{
			dp[i] = ( b(i,0) - (dp[(i-1)]*A(i,(i-1)) ) ) / (A(i,i) - (cp[(i-1)]*A(i,(i-1))) );
		}
	}
	
	//Reverse Sweep
	for (int i=(num_rows-1); i>=0; i--)
	{
		if (i==(num_rows-1))
		{
			this->Data[(i*num_cols)+0] = dp[i];
		}
		else
		{
			this->Data[(i*num_cols)+0] = dp[i] - (cp[i]*this->Data[((i+1)*num_cols)+0]);
		}
	}
	cp.clear();
	dp.clear();
	
	return *this;
}

//Function to solve the tridiagonal matrix using the My Algorithm (works for any tridiagonal matrix)
template <class T>
Matrix<T> &Matrix<T>::ladshawSolve(const Matrix& A, const Matrix& d)
{
	if (this == &A || this == &d)
	{
		mError(arg_matrix_same);
		return *this;
	}
	//Solves the system Ax=b when A is a tridiagonal matrix
	this->set_size(d.num_rows, d.num_cols);
	
	//Check dimensions of A and b to ensure they will work out
	if (num_rows!=A.num_cols || num_cols>1)
	{
		mError(dim_mis_match);
		return *this;
	}
	
	/*
	 * solve for x (i.e. this)
	 * a[i] -> 1 to n-1 (lower diag)
	 * c[i] -> 0 to n-2 (upper diag)
	 * b[i] -> 0 to n-1 (diag)
	 * d[i] -> 0 to n-1 (column vector)
	 * dp = d' ap = a' cp = c' dpp = d" app = a"
	 *
	 * */
	std::vector<T> cp; cp.resize(num_rows);
	std::vector<T> dp; dp.resize(num_rows);
	std::vector<T> ap; ap.resize(num_rows);
	std::vector<T> dpp; dpp.resize(num_rows);
	std::vector<T> app; app.resize(num_rows);
	
	//Forward Sweep
	for (int i=0; i<num_rows; i++)
	{
		if (A(i,i) != 0.0)
		{
			if (i==0)
			{
				ap[i] = 0.0;
				cp[i] = A(i,(i+1))/A(i,i);
			}
			else if (i==(num_rows-1))
			{
				cp[i] = 0.0;
				ap[i] = A(i,(i-1))/A(i,i);
			}
			else
			{
				ap[i] = A(i,(i-1))/A(i,i);
				cp[i] = A(i,(i+1))/A(i,i);
			}
			dp[i] = d(i,0)/A(i,i);
		}
		else
		{
			mError(unstable_matrix);
            mError(singular_matrix);
			ap[i] = 0.0;
			cp[i] = 0.0;
			dp[i] = 0.0;
			return *this;
		}
	}
	
	//Reverse Sweep
	for (int i=(num_rows-1); i>=0; i--)
	{
		if (i==(num_rows-1))
		{
			dpp[i] = dp[i];
			app[i] = ap[i];
		}
		else if (i==0)
		{
			dpp[i] = (dp[i] - (cp[i]*dpp[(i+1)])) / (1 - (cp[i]*app[(i+1)]));
			app[i] = 0;
		}
		else
		{
			dpp[i] = (dp[i] - (cp[i]*dpp[(i+1)])) / (1 - (cp[i]*app[(i+1)]));
			app[i] = ap[i] / (1 - (cp[i]*app[(i+1)]));
		}
	}
	
	//Forward Sweep
	for (int i=0; i<num_rows; i++)
	{
		if (i==0)
			this->Data[(i*num_cols)+0] = dpp[i];
		else
			this->Data[(i*num_cols)+0] = dpp[i] - (app[i] * this->Data[((i-1)*num_cols)+0]);
	}
	cp.clear();
	dp.clear();
	ap.clear();
	dpp.clear();
	app.clear();
	
	return *this;
}


//Function to fill in a tridiagonal matrix with constant coefficients
template <class T>
Matrix<T> &Matrix<T>::tridiagonalFill(const T A, const T B, const T C, bool Spherical)
{
	//Check for square matrix
	if (num_cols!=num_rows)
	{
		mError(non_square_matrix);
		return *this;
	}
	
	//A = lower off diag, B = diag, C = upper off diag
	for (int i=0; i<num_rows; i++)
	{
		for (int j=0; j<num_cols; j++)
		{
			//Diagonals
			if (i==j)
			{
				if (Spherical == false)
					this->Data[(i*num_cols)+j] = B;
				else
					this->Data[(i*num_cols)+j] = B * (j+1);
			}
			//Upper Off Diag
			else if (j==(i+1))
			{
				if (Spherical == false)
					this->Data[(i*num_cols)+j] = C;
				else
					this->Data[(i*num_cols)+j] = C * (j+1);
			}
			//Lower Off Diag
			else if (j==(i-1))
			{
				if (Spherical == false)
					this->Data[(i*num_cols)+j] = A;
				else
					this->Data[(i*num_cols)+j] = A * (j+1);
			}
			else
			{
				this->Data[(i*num_cols)+j] = 0;
			}
		}
	}
	
	return *this;
}

//Fill in a 3D Laplacian matrix with natural ordering of mesh size m
template <class T>
Matrix<T> &Matrix<T>::naturalLaplacian3D(int m)
{
  	int N = m*m*m;
  	this->set_size(N, N);
  	int r = m;
  	int r2 = m*m;
  	int r3 = r2;
  	int r4 = 0;
	
  	for (int i=0; i<N; i++)
  	{
		//If statements for tridiagonal portion
      	this->edit(i, i, 6);
      	if (i == 0)
        {
          	this->edit(i, i+1, -1);
        }
      	else if (i == N-1)
        {
          	this->edit(i, i-1, -1);
        }
      	else if (i == r-1)
        {
          	this->edit(i, i-1, -1);
        }
      	else if (i == r)
        {
          	this->edit(i, i+1, -1);
          	r = r + m;
        }
		else
        {
          	this->edit(i, i-1, -1);
          	this->edit(i, i+1, -1);
        }
		
      	//If statements for 2nd diagonal bands
      	if (i > m-1)
        {
          	if (i <= r3-1)
            {
              	this->edit(i, i-m, -1);
            }
          	else if (i > r3-1)
            {
              	r4 = r4+1;
              	if (r4 == m-1)
                {
                  	r3 = r2;
                  	r4 = 0;
                }
            }
        }
      	if (i <= N-m-1 && i <= r2-m-1)
        {
          	this->edit(i, i+m, -1);
        }
      	if (i == r2-1)
        {
          	r2 = r2+(m*m);
        }
		
      	//If statements for 3rd diagonal bands
      	if (i > (m*m)-1)
        {
          	this->edit(i, i-(m*m), -1);
        }
      	if (i <= N-(m*m)-1)
        {
          	this->edit(i, i+(m*m), -1);
        }
  	}
  	return *this;
}

//Fill in the Dirichlet BC for spherical coordinates
template <class T>
Matrix<T> &Matrix<T>::sphericalBCFill(int node, const T coeff, T variable)
{
	if (num_cols>1)
	{
		mError(dim_mis_match);
		return *this;
	}
	else
	{
		for (int i=0; i<num_rows; i++)
		{
			if (i==(node-1))
			{
				this->Data[(i*num_cols)+0] = (node+1)*coeff*variable;
			}
			else
			{
				this->Data[(i*num_cols)+0] = 0;
			}
		}
		return *this;
	}
}

//Fill in the Constant IC for each system
template <class T>
Matrix<T> &Matrix<T>::ConstantICFill(const T IC)
{
	if (num_cols>1)
	{
		mError(dim_mis_match);
		return *this;
	}
	else
	{
		for (int i=0; i<num_rows; i++)
		{
			this->Data[(i*num_cols)+0] = IC;
		}
		return *this;
	}
}

//Transform Solution Forward or Backward based on flag
template <class T>
Matrix<T> &Matrix<T>::SolnTransform(const Matrix& A, bool Forward)
{
	//Reverse Transform (i.e. q_tilde to q)
	//Forward Transform (i.e. q to q_tilde)
	if (num_cols>1 || num_rows!=A.num_rows)
	{
		mError(dim_mis_match);
		return *this;
	}
	else if (this == &A)
	{
		mError(arg_matrix_same);
		return *this;
	}
	else
	{
		for (int i=0; i<num_rows; i++)
		{
			if (Forward == true)
				this->Data[(i*num_cols)+0] = A.Data[(i*num_cols)+0] * (i+1);
			else
				this->Data[(i*num_cols)+0] = A.Data[(i*num_cols)+0] / (i+1);
		}
		return *this;
	}
}

//Performs Integral average for spherical coordinates based on matrix data (not given center node)
template <class T>
T Matrix<T>::IntegralAvg(double radius, double dr, double bound, bool Dirichlet)
{
	T avg = 0;
	T sum = 0;
	if (num_cols>1)
	{
		mError(dim_mis_match);
		return 0.0;
	}
	for (int i=0; i<num_rows; i++)
	{
		if (i==0)
		{
			//Extrapolation for interior node based on spherical symmetry
			T m, qom, qo;
			m = (this->Data[((i+1)*num_cols)+0] - this->Data[(i*num_cols)+0])/dr;
			qom = this->Data[(i*num_cols)+0] - (m*dr);
			qo = (qom + this->Data[(i*num_cols)+0])/2;
			//Prevent negative mass from bad extrapolation
			if (qo<0)
				qo = 0;
			
			sum = sum + ( ((this->Data[(i*num_cols)+0] + qo )/2)*(pow((i+1),3) - pow(i,3)) );
		}
		else
		{
			sum = sum + ( ((this->Data[(i*num_cols)+0] + this->Data[((i-1)*num_cols)+0])/2)*(pow((i+1),3) - pow(i,3)) );
		}
	}
	//Last sum only needed if using Dirichlet BC at particle edge (if bound is 0, no further sums added)
	if (Dirichlet == true)
		sum = sum + ( ((bound + this->Data[((num_rows-1)*num_cols)+0])/2)*(pow((num_rows+1),3) - pow(num_rows,3)) );
	avg = (pow(dr,3) * sum) / pow(radius,3);
	return avg;
}

//Performs Integral average for spherical coordinates based on matrix data (given center node)
template <class T>
T Matrix<T>::sphericalAvg(double radius, double dr, double bound, bool Dirichlet)
{
	T avg = 0;
	T sum = 0;
	if (num_cols>1)
	{
		mError(dim_mis_match);
		return 0.0;
	}
	for (int i=1; i<num_rows; i++)
	{
		sum = sum + ( ((this->Data[(i*num_cols)+0] + this->Data[((i-1)*num_cols)+0])/2)*(pow((i),3) - pow(i-1,3)) );
	}
	//Last sum only needed if using Dirichlet BC at particle edge (if bound is 0, no further sums added)
	if (Dirichlet == true)
		sum = sum + ( ((bound + this->Data[((num_rows-1)*num_cols)+0])/2)*(pow((num_rows),3) - pow(num_rows-1,3)) );
	avg = (pow(dr,3) * sum) / pow(radius,3);
	return avg;
}

//Performs Integral total for spherical coordinates based on matrix data
template <class T>
T Matrix<T>::IntegralTotal(double dr, double bound, bool Dirichlet)
{
	T tot = 0;
	T sum = 0;
	if (num_cols>1)
	{
		mError(dim_mis_match);
		return 0.0;
	}
	for (int i=0; i<num_rows; i++)
	{
		if (i==0)
		{
			T m, qom, qo;
			m = (this->Data[((i+1)*num_cols)+0] - this->Data[(i*num_cols)+0])/dr;
			qom = this->Data[(i*num_cols)+0] - (m*dr);
			qo = (qom + this->Data[(i*num_cols)+0])/2;
			//Prevent negative mass from bad extrapolation
			if (qo<0)
				qo = 0;
			
			sum = sum + ( ((this->Data[(i*num_cols)+0] + qo )/2)*(pow((i+1),3) - pow(i,3)) );
		}
		else
		{
			sum = sum + ( ((this->Data[(i*num_cols)+0] + this->Data[((i-1)*num_cols)+0])/2)*(pow((i+1),3) - pow(i,3)) );
		}
	}
	//Last sum only needed if using Dirichlet BC at particle edge
	if (Dirichlet == true)
		sum = sum + ( ((bound + this->Data[((num_rows-1)*num_cols)+0])/2)*(pow((num_rows+1),3) - pow(num_rows,3)) );
	tot = (4*M_PI*pow(dr,3)*sum)/3;
	return tot;
}

//Function to fill in the tridiagonal matrix for 1-D Conservation Laws
template <class T>
Matrix<T> &Matrix<T>::tridiagonalVectorFill(const std::vector<T> &A, const std::vector<T> &B, const std::vector<T> &C)
{
	//Check for square matrix
	if (num_cols!=num_rows)
	{
		mError(non_square_matrix);
		return *this;
	}
	if (num_rows!=A.size() || num_rows!=B.size() || num_rows!=C.size())
	{
		mError(matvec_mis_match);
		return *this;
	}
	
  	//A = lower off diag, B = diag, C = upper off diag
	for (int i=0; i<num_rows; i++)
	{
    	//Diagonals
    	this->Data[(i*num_cols)+i] = B[i];
		
    	if (i==0)
    	{
      		//Upper Diagonal
      		this->Data[(i*num_cols)+(i+1)] = C[i];
    	}
    	else if (i==num_rows-1)
    	{
      		//Lower Diagonal
     	 	this->Data[(i*num_cols)+(i-1)] = A[i];
    	}
    	else
    	{
      		//Both
      		this->Data[(i*num_cols)+(i+1)] = C[i];
      		this->Data[(i*num_cols)+(i-1)] = A[i];
    	}
	}
	
	return *this;
}

//Function to fill in a column matrix based on a vector of values
template <class T>
Matrix<T> &Matrix<T>::columnVectorFill(const std::vector<T> &A)
{
	if (num_cols>1)
	{
		mError(dim_mis_match);
		return *this;
	}
	if (num_rows!=A.size())
	{
		mError(matvec_mis_match);
		return *this;
	}
	//Fill in column matrix based on vector of values
	for (int i=0; i<num_rows; i++)
	{
		this->Data[(i*num_cols)+0] = A[i];
	}
	return *this;
}

//Function to perform a column projection based on previous time series info
template <class T>
Matrix<T> &Matrix<T>::columnProjection(const Matrix& b, const Matrix& b_old, const double dt, const double dt_old)
{
	if (b.num_cols>1 || b_old.num_cols>1 || b.num_rows!=b_old.num_rows)
	{
		mError(dim_mis_match);
		return *this;
	}
	if (this == &b || this == &b_old)
	{
		mError(arg_matrix_same);
		return *this;
	}
	this->set_size(b.num_rows,b.num_cols);
	double dtm1 = dt_old;
	if (dtm1 == 0.0) dtm1 = dt;
	for (int i=0; i<num_rows; i++)
	{
		this->Data[(i*num_cols)+0] = b(i,0) + ( (dt/(2.0*dtm1)) * (b(i,0) - b_old(i,0)) );
		//Note: This function does not check for negative mass
	}
	return *this;
}

//Fill in the Dirichlet BC for non-spherical coordinates
template <class T>
Matrix<T> &Matrix<T>::dirichletBCFill(int node, const T coeff, T variable)
{
	if (num_cols>1)
	{
		mError(dim_mis_match);
		return *this;
	}
	else
	{
		for (int i=0; i<num_rows; i++)
		{
			if (i==node)
			{
				this->Data[(i*num_cols)+0] = coeff*variable;
			}
			else
			{
				this->Data[(i*num_cols)+0] = 0;
			}
		}
		return *this;
	}
}

//Directly solve a Diagonal linear system
template <class T>
Matrix<T> &Matrix<T>::diagonalSolve(const Matrix& D, const Matrix& v)
{
  	if (this == &D || this == &v)
  	{
		mError(arg_matrix_same);
		return *this;
  	}
  	//Solves the system Dx=v when D is a diagonal matrix
  	this->set_size(v.num_rows, v.num_cols);
	
  	//Check dimensions of U and v to ensure they will work out
  	if (num_rows!=D.num_cols || num_cols>1)
  	{
		mError(dim_mis_match);
		return *this;
  	}
	
  	//Loop over the diagonals
  	for (int i=0; i<D.num_rows; i++)
  	{
      	if (D(i,i) != 0.0)
        {
      		this->edit(i, 0, v(i,0)/D(i,i));
        }
      	else
        {
          	mError(singular_matrix);
          	return *this;
        }
  	}
  	return *this;
}

//Directly solve an upper triangular linear system
template <class T>
Matrix<T> &Matrix<T>::upperTriangularSolve(const Matrix& U, const Matrix& v)
{
    if (this == &U || this == &v)
    {
        mError(arg_matrix_same);
        return *this;
    }
    //Solves the system Ux=v when U is an upper triangular matrix
    this->set_size(v.num_rows, v.num_cols);
    
    //Check dimensions of U and v to ensure they will work out
    if (num_rows!=U.num_cols || num_cols>1)
    {
        mError(dim_mis_match);
        return *this;
    }
    
    T sum;
    for (int i=U.num_rows-1; i>=0; i--)
    {
        sum = 0.0;
        for (int j=U.num_cols-1; j>i; j--)
        {
            sum = sum + U(i,j) * this->operator()(j, 0);
        }
        if (U(i,i) == 0.0)
        {
            mError(singular_matrix)
            return *this;
        }
        else
        {
            this->edit(i, 0, (v(i,0) - sum) / U(i,i));
        }
    }
    
    return *this;
}

//Directly solve a lower triangular linear system
template <class T>
Matrix<T> &Matrix<T>::lowerTriangularSolve(const Matrix& L, const Matrix& v)
{
    if (this == &L || this == &v)
    {
        mError(arg_matrix_same);
        return *this;
    }
    //Solves the system Ux=v when U is an upper triangular matrix
    this->set_size(v.num_rows, v.num_cols);
    
    //Check dimensions of U and v to ensure they will work out
    if (num_rows!=L.num_cols || num_cols>1)
    {
        mError(dim_mis_match);
        return *this;
    }
    
    T sum;
    for (int i=0; i<L.num_rows; i++)
    {
        sum = 0.0;
        for (int j=0; j<i; j++)
        {
            sum = sum + L(i,j) * this->operator()(j, 0);
        }
        if (L(i,i) == 0.0)
        {
            mError(singular_matrix)
            return *this;
        }
        else
        {
            this->edit(i, 0, (v(i,0) - sum) / L(i,i));
        }
    }
    
    return *this;
}

//Function to convert an upper Hessenberg to upper triangular matrix while editing matrix b
template <class T>
Matrix<T> &Matrix<T>::upperHessenberg2Triangular(Matrix& b)
{
    if (this->num_rows<2 || this->num_cols<2)
    {
        mError(matrix_too_small);
        return *this;
    }
    if (this == &b)
    {
        mError(arg_matrix_same);
        return *this;
    }
    if (this->num_rows != b.num_rows || b.num_cols>1)
    {
        mError(dim_mis_match);
        return *this;
    }
    T s,c;
    T value_i, value_ip1;
    //Loop from top to bottow editing this and b
    for (int i=0; i<this->num_rows-1; i++)
    {
        s = this->operator()(i+1,i) / sqrt( pow(this->operator()(i, i),2.0) + pow(this->operator()(i+1,i),2.0) );
        c = this->operator()(i,i) / sqrt( pow(this->operator()(i,i),2.0) + pow(this->operator()(i+1,i),2.0) );
		
		if (isnan(s) || isnan(c))
		{
			mError(singular_matrix);
			return *this;
		}
        
        value_i = ((b(i,0)*c) + (b(i+1,0)*s));
        value_ip1 = (-(b(i,0)*s) + (b(i+1,0)*c));
        b.edit(i, 0, value_i );
        b.edit(i+1, 0, value_ip1 );
        
        for (int j=i; j<this->num_cols; j++)
        {
            value_i = ((this->operator()(i, j)*c) + (this->operator()(i+1, j)*s));
            value_ip1 = (-(this->operator()(i, j)*s) + (this->operator()(i+1, j)*c));
            this->edit(i, j, value_i);
            this->edit(i+1, j, value_ip1);
        }
    }
    return *this;
}

//Function to convert an lower Hessenberg to lower triangular matrix while editing matrix b
template <class T>
Matrix<T> &Matrix<T>::lowerHessenberg2Triangular(Matrix& b)
{
    if (this->num_rows<2 || this->num_cols<2)
    {
        mError(matrix_too_small);
        return *this;
    }
    if (this == &b)
    {
        mError(arg_matrix_same);
        return *this;
    }
    if (this->num_rows != b.num_rows || b.num_cols>1)
    {
        mError(dim_mis_match);
        return *this;
    }
    T s,c;
    T value_i, value_ip1;
    //Loop from bottom to top editing this and b
    for (int i=this->num_rows-1; i>0; i--)
    {
        s = this->operator()(i-1,i) / sqrt( pow(this->operator()(i, i),2.0) + pow(this->operator()(i-1,i),2.0) );
        c = this->operator()(i,i) / sqrt( pow(this->operator()(i,i),2.0) + pow(this->operator()(i-1,i),2.0) );
		
		if (isnan(s) || isnan(c))
		{
			mError(singular_matrix);
			return *this;
		}
        
        value_i = ((b(i,0)*c) + (b(i-1,0)*s));
        value_ip1 = (-(b(i,0)*s) + (b(i-1,0)*c));
        b.edit(i, 0, value_i );
        b.edit(i-1, 0, value_ip1 );
        
        for (int j=i; j>=0; j--)
        {
            value_i = ((this->operator()(i, j)*c) + (this->operator()(i-1, j)*s));
            value_ip1 = (-(this->operator()(i, j)*s) + (this->operator()(i-1, j)*c));
            this->edit(i, j, value_i);
            this->edit(i-1, j, value_ip1);
        }
    }
    return *this;
}

//Function to solve Hx=v when H is upper Hessenberg
template <class T>
Matrix<T> &Matrix<T>::upperHessenbergSolve(const Matrix& H, const Matrix& v)
{
    if (this == &H || this == &v)
    {
        mError(arg_matrix_same);
        return *this;
    }
    Matrix<T> U(H), b(v);
    U.upperHessenberg2Triangular(b);
    this->upperTriangularSolve(U, b);
    U.Data.clear();
    b.Data.clear();
    return *this;
}

//Function to solve Hx=v when H is lower Hessenberg
template <class T>
Matrix<T> &Matrix<T>::lowerHessenbergSolve(const Matrix& H, const Matrix& v)
{
    if (this == &H || this == &v)
    {
        mError(arg_matrix_same);
        return *this;
    }
    Matrix<T> L(H), b(v);
    L.lowerHessenberg2Triangular(b);
    this->lowerTriangularSolve(L, b);
    L.Data.clear();
    b.Data.clear();
    return *this;
}

//Function to extract a column matrix 'this' from a full matrix M
template <class T>
Matrix<T> &Matrix<T>::columnExtract(int j, const Matrix& M)
{
    if (this == &M)
    {
        mError(arg_matrix_same);
        return *this;
    }
  	if (this->rows() != M.num_rows)
    {
      	this->set_size(M.num_rows, 1);
    }
	
    for (int i=0; i<M.num_rows; i++)
    {
        this->edit(i, 0, M(i,j));
    }
    return *this;
}

//Function to extract a row matrix 'this' from a full matrix M
template <class T>
Matrix<T> &Matrix<T>::rowExtract(int i, const Matrix& M)
{
    if (this == &M)
    {
        mError(arg_matrix_same);
        return *this;
    }
  	if (this->columns() != M.num_cols)
    {
      	this->set_size(1, M.num_cols);
    }
	
    for (int j=0; j<M.num_cols; j++)
    {
        this->edit(0, j, M(i,j));
    }
    return *this;
}

//Function to replace an existing column of 'this' with the column matrix v
template <class T>
Matrix<T> &Matrix<T>::columnReplace(int j, const Matrix &v)
{
    if (this->num_rows != v.num_rows)
    {
        mError(matvec_mis_match)
        return *this;
    }
    if (this == &v)
    {
        mError(arg_matrix_same);
        return *this;
    }
    for (int i=0; i<this->num_rows; i++)
    {
        this->edit(i, j, v(i,0));
    }
    return *this;
}

//Function to replace an existing column of 'this' with the column matrix v
template <class T>
Matrix<T> &Matrix<T>::rowReplace(int i, const Matrix &v)
{
    if (this->num_cols != v.num_cols)
    {
        mError(matvec_mis_match)
        return *this;
    }
    if (this == &v)
    {
        mError(arg_matrix_same);
        return *this;
    }
    for (int j=0; j<this->num_cols; j++)
    {
        this->edit(i, j, v(0,j));
    }
    return *this;
}

//Function to reduce the number of rows in the matrix by neglecting the last row
template <class T>
void Matrix<T>::rowShrink()
{
    /*
     Note: because of how the matrix data is stored, this is the only line
     needed to complete this operation. However, this does not actually remove
     any data. It instead just instructs all other functions to neglect data
     in the last row of the full matrix.
     */
    this->num_rows--;
}

//Function to reduce the number of columns in the matrix by neglecting the last column
template <class T>
void Matrix<T>::columnShrink()
{
    //Here we must force a reinitialization of the data because of how it is stored
    Matrix<T> temp(*this);
    this->set_size(temp.num_rows, temp.num_cols-1);
    for (int i=0; i<this->num_rows; i++)
    {
        for (int j=0; j<this->num_cols; j++)
        {
            this->edit(i, j, temp(i,j));
        }
    }
    temp.Data.clear();
}

//Function to add another row to the end of an existing matrix
template <class T>
void Matrix<T>::rowExtend(const Matrix &v)
{
    if (this->num_cols != v.num_cols)
    {
        mError(matvec_mis_match)
        return;
    }
    for (int j=0; j<v.num_cols; j++)
    {
        this->Data.push_back(v(0,j));
    }
    this->num_rows++;
}

//Function to add another column to the end of an existing matrix
template <class T>
void Matrix<T>::columnExtend(const Matrix &v)
{
    if (this->num_rows != v.num_rows)
    {
        mError(matvec_mis_match)
        return;
    }
    Matrix<T> temp(*this);
    this->set_size(temp.num_rows, temp.num_cols+1);
    for (int i=0; i<this->num_rows; i++)
    {
        for (int j=0; j<this->num_cols; j++)
        {
            if (j<this->num_cols-1)
                this->edit(i, j, temp(i,j));
            else
                this->edit(i, j, v(i,0));
        }
    }
}

int MACAW_TESTS();


#endif /* MACAW_HPP_ */
