#include"utils.h"
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<float.h>

const float EPS = 1e-6;

void initMatrixWithRV(float *A, int rows, int cols)
{
    srand((unsigned)time(NULL));
    for(int i = 0; i < rows*cols; i++){
        A[i] = (float)rand() / RAND_MAX;
    }
}

void initMatrixWithZero(float *A, int rows, int cols)
{
    for(int i = 0; i < rows*cols; i++){
        A[i] = 0;
    }
}

void printMatrix(float *A, int rows, int cols)
{
    if(cols > 10){
        std::cout << "Too many columns. Can not show properly. " << std::endl;
        return;
    }
    
    std::cout.precision(2);
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            std::cout << A[i*cols + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void matrix_multiply(float *A, float *B, float *matResult, int m, int p, int n)
{
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            float temp = 0;
            for(int k=0; k<p; k++){
                temp += A[i*p+k] * B[k*n + j];
            }
            matResult[i*n+j] = temp;
        }
    }
}

void matrix_multiply_with_tB(float *A, float *tB, float *matResult, int m, int n, int p)
{
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            float temp = 0;
            for(int k=0; k<p; k++){
                temp += A[i*p+k] * tB[j*p+k];
            }
            matResult[i*n+j] = temp;
        }
    }
}

bool isMatrixEqual(float *A, float *B, int rows, int cols)
{
    const float EPS = 1e-6;
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            if(fabs(A[i*cols+j] - B[i*cols+j]) > EPS){
                std::cout << "The id is : i = " << i << " j = " << j << std::endl;
                std::cout << "The gap is : " << fabs(A[i*cols+j] - B[i*cols+j]) << std::endl;
                return false;
            }
        }
    }
    
    return true;
}

float max(float *vec, int length)
{
    float maximum = FLT_MIN;
    for(int i=0; i<length; i++){
        if(vec[i] > maximum){
            maximum = vec[i];
        }
    }
    
    return maximum;
}

void copyMatrix(float *A, float *A_copy, int rows, int cols)
{
    for(int i=0; i<rows*cols; i++){
        A_copy[i] = A[i];
    }
}

void transpose_matrix(float *A, int rows, int cols){
    float* B = new float[rows*cols];
    for(int i=0; i<rows*cols; i++){
        B[i] = A[i];
    }
    
    for(int i=0; i<cols; i++){
        for(int j=0; j<rows; j++){
            A[i*rows+j] = B[j*cols+i];
        }
    }
    delete[] B;
}





