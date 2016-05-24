#include"SCVector.h"
#include<iostream>
#include<math.h>

SCVector::SCVector(int dim){
	dimension = dim;
	data = new double[dim];
	for(int i=0; i<dimension; i++){
		data[i] = 0.0;
	}
}

SCVector::SCVector(const SCVector& v){
	dimension = v.dimension;
	data = new double[dimension];
	for(int i=0; i<dimension; i++){
		data[i] = v(i);
	}
}

SCVector::SCVector(int col, const SCMatrix& A){
	dimension = A.getRows();
	data = new double[dimension];
	for(int i=0; i<A.getRows(); i++){
		data[i] = A(i, col);
	}
}

SCVector::~SCVector(){
	delete[] data;
	data = NULL;
	dimension = 0;
}

int SCVector::Dimension() const{
	return dimension;
}

double SCVector::Length(){
	double sum = 0;
	for(int i=0; i<dimension; i++){
		sum += data[i]*data[i];
	}
	sum = sqrt(sum);

	return sum;
}

void SCVector::Normalize(){
	double len = Length();
	if(len<1e-6){
		return;
	}

	for(int i=0; i<dimension; i++){
		data[i] /= len;
	}
}

double SCVector::Norm_l1(){
	double sum = 0;
	for(int i=0; i<dimension; i++){
		sum += fabs(data[i]);
	}

	return sum;
}

double SCVector::Norm_l2(){
	return SCVector::Length();
}

double SCVector::Norm_linf(){
	double result = 0;
	for(int i=0; i<dimension; i++){
		if(fabs(data[i]) > result){
			result = fabs(data[i]);
		}
	}

	return result;
}

int SCVector::operator==(const SCVector& v) const{
	if(dimension  != v.Dimension()){
		return 0;
	}

	for(int i=0; i<dimension; i++){
		if(fabs(data[i]-v(i)) > 1e-6){
			return 0;
		}
	}

	return 1;
}

int SCVector::operator!=(const SCVector& v) const{
	return 1-(*this==v);
}

SCVector& SCVector::operator=(const SCVector& v){
	delete[] data;
	data = NULL;	
	dimension = v.Dimension();
	data = new double[dimension];
	for(int i=0; i<dimension; i++){
		data[i] = v(i);
	}

	return *this;
}

double SCVector::operator()(const int i) const{
	return data[i];
}

double& SCVector::operator()(const int i){
	return data[i];
}

void SCVector::Print() const{
	for(int i=0; i<dimension; i++){
		std::cout << data[i] << " ";
		if(i%10 == 9){
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}

void SCVector::Initialize(double a){
}

void SCVector::Initialize(double* v){
}


