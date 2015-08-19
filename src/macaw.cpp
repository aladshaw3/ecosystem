//----------------------------------------
//  Created by Austin Ladshaw on 1/7/14
//  Copyright (c) 2014
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*
 *		MACAW = MAtrix CAlculation Workshop
 *
 *		This is a small C++ library that faciltates the use and construction of
 *		real matrices using vector objects. It has functions for both dense and
 *		sparse matrices, though the sparse system is poorly supported at current.
 */

#include "macaw.h"

Matrix::Matrix(int rows, int columns)
:
num_rows(rows),
num_cols(columns),
Data(rows * columns)
{

}

//For setting values
double& Matrix::operator()(int i, int j)
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
double Matrix::operator()(int i, int j) const
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
Matrix::Matrix(const Matrix& M)
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
Matrix &Matrix::operator=(const Matrix &M)
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
Matrix::Matrix()
:
Data(1)
{
	num_cols = 1;
	num_rows = 1;
}

//Default Destructor
Matrix::~Matrix()
{
    Data.clear();
}

//Allocate size
void Matrix::set_size(int i, int j)
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
void Matrix::zeros()
{
  	for (int n=0; n<this->Data.size(); n++)
    		this->Data[n] = 0.0;
}

//Add/Change an element to a sparse matrix
void Matrix::edit(int i, int j, double value)
{
    if (i>=num_rows || j>=num_cols || i<0 || j<0)
    {
        mError(out_of_bounds);
        return;
    }
    this->Matrix::operator()(i, j) = value;
}

//Function to return the number of rows in the matrix
int Matrix::rows()
{
    return this->num_rows;
}

//Function to return the number of columns in the matrix
int Matrix::columns()
{
    return this->num_cols;
}

//Function to determine the determinate of a matrix
double Matrix::determinate()
{
	double det=0;
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
			Matrix temp((num_rows-1), (num_cols-1));
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
double Matrix::norm()
{
  double norm = 0.0;
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
double Matrix::sum()
{
	double sum = 0.0;
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
double Matrix::inner_product(const Matrix& x)
{
  double prod = 0.0;
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
Matrix &Matrix::cofactor(const Matrix &M)
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
			Matrix temp((M.num_rows-1), (M.num_cols-1));

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
Matrix Matrix::operator+(const Matrix &M)
{
	if(num_rows!=M.num_rows && num_cols!=M.num_cols)
	{
		mError(dim_mis_match);
		return *this;
	}
	else
	{
		Matrix temp(num_rows,num_cols);
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
Matrix Matrix::operator-(const Matrix &M)
{
	if(num_rows!=M.num_rows && num_cols!=M.num_cols)
	{
		mError(dim_mis_match);
		return *this;
	}
	else
	{
		Matrix temp(num_rows,num_cols);
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
Matrix Matrix::operator*(const double a)
{
	Matrix temp(num_rows,num_cols);
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
Matrix Matrix::operator/(const double a)
{
	Matrix temp(num_rows,num_cols);
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
Matrix Matrix::operator*(const Matrix &M)
{
	if (num_cols!=M.num_rows)
	{
		mError(dim_mis_match);
		return *this;
	}
	else
	{
		Matrix temp(num_rows, M.num_cols);
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
Matrix &Matrix::transpose(const Matrix &M)
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
Matrix &Matrix::transpose_multiply(const Matrix &MT, const Matrix &v)
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
Matrix &Matrix::adjoint(const Matrix &M)
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
Matrix &Matrix::inverse(const Matrix &M)
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
		Matrix A(M);
		det_inv = 1/A.determinate();
		*this = this->adjoint(A) * det_inv;
		A.Data.clear();
		return *this;
	}
}

//Display Matrix to console
void Matrix::Display(const std::string Name)
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
Matrix &Matrix::tridiagonalSolve(const Matrix& A, const Matrix& b)
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
	std::vector<double> cp; cp.resize(num_rows-1);
	std::vector<double> dp; dp.resize(num_rows);

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
Matrix &Matrix::ladshawSolve(const Matrix& A, const Matrix& d)
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
	std::vector<double> cp; cp.resize(num_rows);
	std::vector<double> dp; dp.resize(num_rows);
	std::vector<double> ap; ap.resize(num_rows);
	std::vector<double> dpp; dpp.resize(num_rows);
	std::vector<double> app; app.resize(num_rows);

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
Matrix &Matrix::tridiagonalFill(const double A, const double B, const double C, bool Spherical)
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
Matrix &Matrix::naturalLaplacian3D(int m)
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
Matrix &Matrix::sphericalBCFill(int node, const double coeff, double variable)
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
Matrix &Matrix::ConstantICFill(const double IC)
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
Matrix &Matrix::SolnTransform(const Matrix& A, bool Forward)
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
double Matrix::IntegralAvg(double radius, double dr, double bound, bool Dirichlet)
{
	double avg = 0;
	double sum = 0;
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
			double m, qom, qo;
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
double Matrix::sphericalAvg(double radius, double dr, double bound, bool Dirichlet)
{
	double avg = 0;
	double sum = 0;
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
double Matrix::IntegralTotal(double dr, double bound, bool Dirichlet)
{
	double tot = 0;
	double sum = 0;
	if (num_cols>1)
	{
		mError(dim_mis_match);
		return 0.0;
	}
	for (int i=0; i<num_rows; i++)
	{
		if (i==0)
		{
			double m, qom, qo;
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
Matrix &Matrix::tridiagonalVectorFill(const std::vector<double> &A, const std::vector<double> &B, const std::vector<double> &C)
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
Matrix &Matrix::columnVectorFill(const std::vector<double> &A)
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

Matrix &Matrix::columnProjection(const Matrix& b, const Matrix& b_old, const double dt, const double dt_old)
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
Matrix &Matrix::dirichletBCFill(int node, const double coeff, double variable)
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
Matrix &Matrix::diagonalSolve(const Matrix& D, const Matrix& v)
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
Matrix &Matrix::upperTriangularSolve(const Matrix& U, const Matrix& v)
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
    
    double sum;
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
Matrix &Matrix::lowerTriangularSolve(const Matrix& L, const Matrix& v)
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
    
    double sum;
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
Matrix &Matrix::upperHessenberg2Triangular(Matrix& b)
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
    double s,c;
    double value_i, value_ip1;
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
Matrix &Matrix::lowerHessenberg2Triangular(Matrix& b)
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
    double s,c;
    double value_i, value_ip1;
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
Matrix &Matrix::upperHessenbergSolve(const Matrix& H, const Matrix& v)
{
    if (this == &H || this == &v)
    {
        mError(arg_matrix_same);
        return *this;
    }
    Matrix U(H), b(v);
    U.upperHessenberg2Triangular(b);
    this->upperTriangularSolve(U, b);
    U.Data.clear();
    b.Data.clear();
    return *this;
}

//Function to solve Hx=v when H is lower Hessenberg
Matrix &Matrix::lowerHessenbergSolve(const Matrix& H, const Matrix& v)
{
    if (this == &H || this == &v)
    {
        mError(arg_matrix_same);
        return *this;
    }
    Matrix L(H), b(v);
    L.lowerHessenberg2Triangular(b);
    this->lowerTriangularSolve(L, b);
    L.Data.clear();
    b.Data.clear();
    return *this;
}

//Function to extract a column matrix 'this' from a full matrix M
Matrix &Matrix::columnExtract(int j, const Matrix& M)
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
Matrix &Matrix::rowExtract(int i, const Matrix& M)
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
Matrix &Matrix::columnReplace(int j, const Matrix &v)
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
Matrix &Matrix::rowReplace(int i, const Matrix &v)
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
void Matrix::rowShrink()
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
void Matrix::columnShrink()
{
    //Here we must force a reinitialization of the data because of how it is stored
    Matrix temp(*this);
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
void Matrix::rowExtend(const Matrix &v)
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
void Matrix::columnExtend(const Matrix &v)
{
    if (this->num_rows != v.num_rows)
    {
        mError(matvec_mis_match)
        return;
    }
    Matrix temp(*this);
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

int MACAW_TESTS()
{
    int success = 0;
    
    //Test Lower Triangular solve (PASS)
    int rows, cols;
    rows = 5;
    cols = 5;
    Matrix L(rows,cols);
    Matrix v(rows,1);
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<=i; j++)
        {
            L.edit(i,j, i+j+1);
        }
        v.edit(i, 0, i+1);
    }
    L.Display("L");
    v.Display("v");
    Matrix x;
    x.lowerTriangularSolve(L, v);
    x.Display("x");
    (L*x).Display("L*x=v=");
    
    //Test of Upper Triangular solve (PASS)
    Matrix U(rows,cols);
    for (int i=rows-1; i>=0; i--)
    {
        for (int j=cols-1; j>=i; j--)
        {
            U.edit(i,j, i+j+1);
        }
    }
    U.Display("U");
    x.upperTriangularSolve(U, v);
    x.Display("x");
    (U*x).Display("U*x=v=");
    
    //Test row and column reductions
    U.rowShrink(); //(PASS)
    U.Display("U");
    U.columnShrink(); //(PASS - had to reinitialize matrix to work)
    U.Display("U");
    
    //Test Hessenberg Functions (PASS)
    Matrix H(rows,cols);
    for (int i=rows-1; i>=0; i--)
    {
        for (int j=cols-1; j>=i; j--)
        {
            H.edit(i,j, i+j+1);
            if (i>=0 && i<rows-1)
            {
                //H.edit(i+1,j, i+j+1); //Matrix is singular!!!
                H.edit(i+1,j, i+i+1); //Note: changed j to i to prevent singularity
            }
        }
    }
    H.Display("H");
    v.Display("v");
    Matrix H2, v2;
    H2 = H;
    v2 = v;
    H2.upperHessenberg2Triangular(v2);
    H2.Display("H2");
    v2.Display("v2");
    x.upperTriangularSolve(H2, v2);
    x.Display("x");
    (H*x).Display("should be v");
    
    Matrix H3(rows,cols);
    H3.tridiagonalFill(-1, 2, -1, true);
    H3.Display("H3");
    Matrix H4, v4;
    H4 = H3;
    v4 = v;
    H4.lowerHessenberg2Triangular(v4);
    H4.Display("H4");
    v4.Display("v4");
    x.lowerTriangularSolve(H4, v4);
    x.Display("x");
    x.ladshawSolve(H3, v);
    x.Display("x");
    
    x.upperHessenbergSolve(H, v);
    x.Display("x");
    (H*x).Display("v");
    
    x.lowerHessenbergSolve(H3, v);
    x.Display("x");
    (H3*x).Display("v");
    
    H.Display("H");
    H3.Display("H3");
    
    //Warning! Conversion of a Hessenberg to Triangular matrix may result in singularity!!!
    
    //Test the row and column extract and replace functions (PASS)
    Matrix r, c;
    r.rowExtract(2, H3);
    r.Display("r");
    c.columnExtract(2, H3);
    c.Display("c");
    Matrix H5;
    H5 = H3;
    H5.columnReplace(0, c);
    H5.Display("H5");
    H5.rowReplace(0, r);
    H5.Display("H5");
    H5.rowReplace(0, r.transpose(c));
    H5.Display("H5");
    
    //Test the column and row extension functions (PASS)
    Matrix A(rows,cols);
    A.tridiagonalFill(-1, 2, -1, false);
    A.Display("A");
    Matrix cex(rows,1);
    Matrix rex(1,cols);
    cex.dirichletBCFill(rows-1, 2, 1);
    cex.Display("cex");
    rex.transpose(cex);
    rex.Display("rex");
    
    cex.set_size(rows+1, 1);
    cex.dirichletBCFill(rows, 1, 2);
    A.rowExtend(rex);
    A.Display("A");
    
    //Test of Upper Hessenberg non-square conversion (PASS)
    Matrix HA;
    HA = A;
    HA.upperHessenberg2Triangular(cex);
    HA.Display("HA");
    
    HA.rowShrink();
    HA.Display("HA_Tri");
    
    A.columnExtend(cex);
    A.Display("A");
  
  	//Test of Diagonal Solver
  	Matrix D(10,10), b(10,1);
  	D.tridiagonalFill(0, 2, 0, true);
  	D.Display("D");
  	b.ConstantICFill(1.0);
  	Matrix xD;
  	xD.diagonalSolve(D, b);
  	xD.Display("x");
  	std::cout << "norm = " << (b - D*xD).norm() << std::endl;
  
  
    return success;
}


