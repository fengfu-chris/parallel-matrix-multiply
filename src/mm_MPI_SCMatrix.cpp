/* parallel matrix multiplication with MPI and SCMatrix
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
#include"SCMatrix.h"

const double EPS = 1e-6;
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

	
    double elapseTime, elapseTime_single; 

    int myrank, numprocs;

    MPI_Status status;
  
	
    MPI_Init(&argc, &argv);  // 并行开始
	//std::cout << "Start to computing..." << std::endl;
  
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  // 获取进程数
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  // 获取本进城ID号
 
    if(numprocs < 2){
        std::cout << "Error! There must be more than two processors." << std::endl;
		return 1;
    }
    
    int bm = m / numprocs;
    int bn = n / numprocs;
     
	SCMatrix A(1), B(1), C(1), C_true(1);
	SCMatrix bA(bm,p), bB_send(bn,p), bB_recv(bn,p), bC(bm, bn), bC_send(bm, n);
 	  
    if(myrank == 0){
		SCMatrix A1(m,p), B1(n,p), C1(m,n), C_true1(m,n);

		A = A1; B = B1; C = C1; C_true = C_true1;

		A.init_with_RV();
		B.init_with_RV();    
    }
    MPI_Barrier(MPI_COMM_WORLD);
	//A.print(); B.print();
    
	
    // start to timing
    clock_t start_time = clock();    

    //clock_t t_temp_start = clock();
    MPI_Scatter(&A, bm * p, MPI_DOUBLE, &bA, bm * p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(&B, bn * p, MPI_DOUBLE, &bB_recv, bn * p, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int sendTo = (myrank + 1) % numprocs;
	int recvFrom = (myrank - 1 + numprocs) % numprocs;
		
	int circle = 0;
	
	bB_send = bB_recv;
	do{
		bC = bA * bB_recv.transpose();
		
		int blocks_col = (myrank - circle + numprocs) % numprocs;
		for(int i=0; i<bm; i++){
			for(int j=0; j<bn; j++){
				bC_send(i, blocks_col*bn + j) = bC(i, j);
			}
		}

		if(myrank % 2 == 0){
			bB_send = bB_recv;
			MPI_Ssend(&bB_send, bn*p, MPI_DOUBLE, sendTo, circle, MPI_COMM_WORLD);
			MPI_Recv(&bB_recv, bn*p, MPI_DOUBLE, recvFrom, circle, MPI_COMM_WORLD, &status);
		}else{
			MPI_Recv(&bB_recv, bn*p, MPI_DOUBLE, recvFrom, circle, MPI_COMM_WORLD, &status);
			MPI_Ssend(&bB_send, bn*p, MPI_DOUBLE, sendTo, circle, MPI_COMM_WORLD);
			bB_send = bB_recv;
		}
		
		circle++;
    }while(circle < numprocs);

	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Gather(&bC_send, bm * n, MPI_DOUBLE, &C, bm * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
    if(myrank == 0){
		int remainAStartId = bm * numprocs;
		int remainBStartId = bn * numprocs;
		
		for(int i=remainAStartId; i<m; i++){
			for(int j=0; j<n; j++){
				double temp=0;
				for(int k=0; k<p; k++){
					temp += A(i,k) * B(j,k);
				}
				C(i,j) = temp;
			}
		}
		
		for(int i=0; i<remainAStartId; i++){
			for(int j=remainBStartId; j<n; j++){
				double temp = 0;
				for(int k=0; k<p; k++){
					temp += A(i,k)*B(j,k);
				}
				C(i,j) = temp;
			}
		}
    }
    
    // end timing
    clock_t end_time = clock();	

    elapseTime = (double)(end_time-start_time) / CLOCKS_PER_SEC;  

    if(myrank == 0){
		//std::cout << "[P_0] m = " << m << ", np = " << numprocs << ", time cost: " << elapseTime  << std::endl;
	    std::cout << "[P_0] Totally cost " << elapseTime << " seconds in parallel procedure." << std::endl;
        
		start_time = clock();
		C_true = A * B.transpose();
        end_time = clock();

		elapseTime_single = double(end_time-start_time) / CLOCKS_PER_SEC;
		std::cout << "[P_0] Totally cost " << elapseTime_single << " seconds with one processor." << std::endl;

        bool b = (C==C_true);
        if(b){
            std::cout << "[P_0] Congradulations! The two results are equal." << std::endl;
        }else{
            std::cout << "[P_0] Bad news! The two results do not match." << std::endl;
			A.print();
			B.print();
			C.print();
			C_true.print();
       }
    }
	
     
    MPI_Finalize(); // 并行结束

    return 0;
}
