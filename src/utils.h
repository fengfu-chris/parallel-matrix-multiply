void initMatrixWithRV(float *A, int rows, int cols);
void initMatrixWithZero(float *A, int rows, int cols);
void printMatrix(float *A, int rows, int cols);
void matrix_multiply(float *A, float *B, float *matResult, int m, int p, int n);
void matrix_multiply_with_tB(float *A, float *tB, float *matResult, int m, int n, int p);
bool isMatrixEqual(float *A, float *B, int rows, int cols);
float max(float *vec, int length);
void copyMatrix(float *A, float *A_copy, int rows, int cols);
void transpose_matrix(float *A, int rows, int cols);
