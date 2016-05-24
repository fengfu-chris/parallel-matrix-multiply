/* parallel matrix multiplication with MPI: updated
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
        std::cout << "Error! Three arguments m, n and p are needed!" << std::endl;
		return 1;
    }
	
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int p = atoi(argv[3]);    

    float *A, *B, *C, *C_true;
    float *bA, *bB_send, *bB_recv, *bC, *bC_send;
    float elapseTime, elapseTime_single;

    int myrank, numprocs;

    MPI_Status status;
  
    MPI_Init(&argc, &argv);  // 并行开始
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  // 获取进程数
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  // 获取本进城ID号
 
    if(numprocs < 2){
        std::cout << "Error! There must be more than two processors." << std::endl;
		return 1;
    }
    
    int bm = m / numprocs;
    int bn = n / numprocs;

    bA = new float[bm * p];
    bB_send = new float[bn * p];
	bB_recv = new float[bn * p];
    bC = new float[bm * bn];
	bC_send = new float[bm * n];
    initMatrixWithZero(bA, bm, p);
    initMatrixWithZero(bB_send, bn, p);
	initMatrixWithZero(bB_recv, bn, p);
    initMatrixWithZero(bC, bm, bn);
	initMatrixWithZero(bC_send, bm, n);
    
    if(myrank == 0){
		//std::cout << "[P_0] m = " << m << ", np = " << numprocs << std::endl;
        A = new float[m * p];
		B = new float[n * p];
        C = new float[m * n];
        C_true = new float[m * n];
        
        initMatrixWithRV(A, m, p);
        initMatrixWithRV(B, n, p);
        initMatrixWithZero(C, m, n);
        initMatrixWithZero(C_true, m, n);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // start to timing
    clock_t start_time = clock();    

    //clock_t t_temp_start = clock();
    MPI_Scatter(A, bm * p, MPI_FLOAT, bA, bm * p, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, bn * p, MPI_FLOAT, bB_recv, bn * p, MPI_FLOAT, 0, MPI_COMM_WORLD);

	int sendTo = (myrank + 1) % numprocs;
	int recvFrom = (myrank - 1 + numprocs) % numprocs;
		
	int circle = 0;
	
	copyMatrix(bB_recv, bB_send, bn, p);
	do{
		matrix_multiply_with_tB(bA, bB_recv, bC, bm, bn, p);
		int blocks_col = (myrank - circle + numprocs) % numprocs;
		for(int i=0; i<bm; i++){
			for(int j=0; j<bn; j++){
				bC_send[i*n + blocks_col*bn + j] = bC[i*bn + j];
			}
		}

		if(myrank % 2 == 0){
			copyMatrix(bB_recv, bB_send, bn, p);
			MPI_Ssend(bB_send, bn*p, MPI_FLOAT, sendTo, circle, MPI_COMM_WORLD);
			//std::cout << "[P_" << myrank << "] step " << circle << " send message to : " << sendTo << std::endl;
			MPI_Recv(bB_recv, bn*p, MPI_FLOAT, recvFrom, circle, MPI_COMM_WORLD, &status);
			//std::cout << "[P_" << myrank << "] step " << circle << " received message from: " << recvFrom << std::endl;
		}else{
			MPI_Recv(bB_recv, bn*p, MPI_FLOAT, recvFrom, circle, MPI_COMM_WORLD, &status);
			//std::cout << "[P_" << myrank << "] step " << circle << " received message from: " << recvFrom << std::endl;	
			MPI_Ssend(bB_send, bn*p, MPI_FLOAT, sendTo, circle, MPI_COMM_WORLD);
			//std::cout << "[P_" << myrank << "] step " << circle << " send message to : " << sendTo << std::endl;	
			copyMatrix(bB_recv, bB_send, bn, p);	
		}
		
		circle++;
    }while(circle < numprocs);

	//std::cout << "[P_" << myrank << "] reached here." << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Gather(bC_send, bm * n, MPI_FLOAT, C, bm * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
    if(myrank == 0){
		int remainAStartId = bm * numprocs;
		int remainBStartId = bn * numprocs;
		
		for(int i=remainAStartId; i<m; i++){
			for(int j=0; j<n; j++){
				float temp=0;
				for(int k=0; k<p; k++){
					temp += A[i*p + k] * B[j*p +k];
				}
				C[i*p + j] = temp;
			}
		}
		
		for(int i=0; i<remainAStartId; i++){
			for(int j=remainBStartId; j<n; j++){
				float temp = 0;
				for(int k=0; k<p; k++){
					temp += A[i*p + k] * B[j*p +k];
				}
				C[i*p + j] = temp;
			}
		}
    }
    
    // end timing
    clock_t end_time = clock();	

    elapseTime = (float)(end_time-start_time) / CLOCKS_PER_SEC;  

    if(myrank == 0){
		//std::cout << "[P_0] m = " << m << ", np = " << numprocs << ", time cost: " << elapseTime  << std::endl;
		std::cout << "[P_0] Totally cost " << elapseTime << " seconds in parallel procedure." << std::endl;

        start_time = clock();
        matrix_multiply_with_tB(A, B, C_true, m, n, p);
        end_time = clock();

		elapseTime_single = float(end_time-start_time) / CLOCKS_PER_SEC;
        
        std::cout << "[P_0] Totally cost " << elapseTime_single << " seconds with one processor." << std::endl;

        bool b = isMatrixEqual(C, C_true, m, n);
        if(b){
            std::cout << "[P_0] Congradulations! The two results are equal." << std::endl;
        }else{
            std::cout << "[P_0] Bad news! The two results do not match." << std::endl;
			//printMatrix(C, m, n);
			//printMatrix(C_true, m, n);
       }
    }
     
    delete[] bA;
    delete[] bB_send;
	delete[] bB_recv;
    delete[] bC;
	delete[] bC_send;
    
    if(myrank == 0){
        delete[] A;
		delete[] B;
        delete[] C;
        delete[] C_true;
    }
    
    MPI_Finalize(); // 并行结束

    return 0;
}
