/* demo of parallel matrix multiplication with MPI
 * C(m,n) = A(m,p) * B(p,n)
 * input: three parameters - m, p, n
 * @copyright: fengfu-chris 2015-6-1
 */

#include<iostream>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>
#include<float.h>
#include"utils.h"

const float EPS = 1e-6;
const int MAX_SPEEDUP = INT_MAX;

int main(int argc, char** argv)
{
    if(argc != 4){
        std::cout << "Error! Three arguments m, p and n are needed!" << std::endl;
    }
    int m = atoi(argv[1]);
    int p = atoi(argv[2]);
    int n = atoi(argv[3]);
    
    float *A, *B, *C, *C_true;
    float *bA, *bC;
    float elapseTime, elapseTimeWithSingleProcessor, *elapseTimeRecv;
 
    int myrank, numprocs;

    MPI_Status status;
  
    MPI_Init(&argc, &argv);  // 并行开始
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  // 获取进程数
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  // 获取本进城ID号

    //std::cout << "[P_" << myrank << "] m = " << m << ", p = " << p << ", n = " << n << std::endl;
    
    if(numprocs < 2){
        std::cout << "Error! There must be more than two processors." << std::endl;
    }
    
    int bm = m / numprocs;

    bA = new float[bm * p];
    B  = new float[p * n];
    bC = new float[bm * n];
    initMatrixWithZero(bA, bm, p);
    initMatrixWithZero(B, p, n);
    initMatrixWithZero(bC, bm, n);
    
    if(myrank == 0){
	elapseTimeRecv = new float[numprocs];
        A = new float[m * p];
        C = new float[m * n];
        C_true = new float[m * n];
        
        initMatrixWithRV(A, m, p);
        initMatrixWithRV(B, p, n);
        initMatrixWithZero(C, m, n);
        initMatrixWithZero(C_true, m, n);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // start to timing
    clock_t start_time = clock();    

    clock_t t_temp_start = clock();
    MPI_Scatter(A, bm * p, MPI_FLOAT, bA, bm *p, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, p * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    clock_t t_temp_end = clock();

    matrix_multiply(bA, B, bC, bm, p, n);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Gather(bC, bm * n, MPI_FLOAT, C, bm * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
    int remainRowsStartId = bm * numprocs;
    if(myrank == 0 && remainRowsStartId < m){
        int remainRows = m - remainRowsStartId;
        matrix_multiply(A + remainRowsStartId * p, B, C + remainRowsStartId * n, remainRows, p, n);
    }
    
    // end timing
    clock_t end_time = clock();	

    elapseTime = (float)(end_time-start_time) / CLOCKS_PER_SEC;
        
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&elapseTime, 1, MPI_FLOAT, elapseTimeRecv, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if(myrank == 0){
        elapseTime = max(elapseTimeRecv, numprocs);
		std::cout << "[P_0] Totally cost " << elapseTime << " seconds in parallel procedure." << std::endl;

        start_time = clock();
        matrix_multiply(A, B, C_true, m, p, n);
        end_time = clock();

		elapseTimeWithSingleProcessor = float(end_time-start_time) / CLOCKS_PER_SEC;
		std::cout << "[P_0] Totally cost " << elapseTimeWithSingleProcessor << " seconds with single processor." << std::endl;

        bool b = isMatrixEqual(C, C_true, m, n);
        if(b){
            std::cout << "[P_0] Congradulations! The two results are equal." << std::endl;
        }else{
            std::cout << "[P_0] Bad news! The two results do not match." << std::endl;

        }
    }
     
    delete[] bA;
    delete[] B;
    delete[] bC;
    
    if(myrank == 0){
		delete[] elapseTimeRecv;
        delete[] A;
        delete[] C;
        delete[] C_true;
    }
    
    MPI_Finalize(); // 并行结束

    return 0;
}   







