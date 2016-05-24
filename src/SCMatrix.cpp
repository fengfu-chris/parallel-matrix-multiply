#include"SCMatrix.h"
#include<iostream>
#include<math.h>
#include<stdlib.h>

SCMatrix::SCMatrix(){
}

SCMatrix::SCMatrix(int m){
	rows = m; 
	cols = m;
	data = new double[rows*cols];
	SCMatrix::init_with_Zero();
}

SCMatrix::SCMatrix(int m, int n){
	rows = m;
	cols = n;
	data = new double[rows*cols];
	SCMatrix::init_with_Zero();
}

SCMatrix::SCMatrix(const SCMatrix& A){
	rows = A.getRows();
	cols = A.getCols();
	data = new double[rows*cols];
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			data[i*cols+j] = A(i,j);
		}
	}
}

SCMatrix::~SCMatrix(){
	delete[] data;
	data = NULL;
	rows = 0;
	cols = 0;
}

int SCMatrix::getRows() const{
	return rows;
}

int SCMatrix::getCols() const{
	return cols;
}

int SCMatrix::operator==(const SCMatrix& A) const{
	if(rows != A.getRows() || cols != A.getCols()){
		return 0;
	}

	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			if(fabs(data[i*cols+j] - A(i,j))>1e-6){
				return 0;
			}
		}
	}
	
	return 1;
}

int SCMatrix::operator!=(const SCMatrix& A) const{
	return 1-(*this==A);
}

SCMatrix& SCMatrix::operator=(const SCMatrix& A){
	delete[] data;
	data = NULL;
	rows = A.getRows();
	cols = A.getCols();
	data = new double[rows*cols];
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			data[i*cols+j] = A(i,j);
		}
	}

	return *this;
}

double* SCMatrix::operator&() const{
	return data;
}

SCMatrix SCMatrix::operator*(const SCMatrix& A){
	SCMatrix B(rows, A.getCols());
	for(int i=0; i<rows; i++){
		for(int j=0; j<A.getCols(); j++){
			double temp = 0;
			for(int k=0; k<cols; k++){
				temp += (*this)(i,k) * A(k,j);
			}
			B(i,j) = temp;
		}
	}

	return B;
}

SCMatrix SCMatrix::transpose() const{
	SCMatrix A(cols, rows);
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			A(j,i) = (*this)(i,j);
		}
	}

	return A;
}

double SCMatrix::operator()(const int i, const int j) const{
	return data[i*cols + j];
}

double& SCMatrix::operator()(const int i, const int j){
	return data[i*cols + j];
}

void SCMatrix::print() const{
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			std::cout << data[i*cols+j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void SCMatrix::init_with_RV(){
    srand((unsigned)time(NULL));
    for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			data[i*cols+j] = (float)rand() / RAND_MAX;
		}
    }
}

void SCMatrix::init_with_Zero(){
    for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			data[i*cols+j] = 0.0;
		}
    }
}
