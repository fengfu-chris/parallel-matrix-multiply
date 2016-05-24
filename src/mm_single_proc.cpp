/* demo of matrix multiplication with single processor
 * C(m,n) = A(m,p) * B(p,n)
 * input: three parameters - m, p, n, repeats
 * @copyright: fengfu-chris 2015-6-1
 */

#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"utils.h"

int main(int argc, char** argv)
{
	std::cout << "Please input three arguments m,p,n and repeat times (optional)." << std::endl;
    if(argc < 4){
        std::cout << "Error! At least three arguments m, p and n are needed!" << std::endl;
    }
    
    int m = atoi(argv[1]);
    int p = atoi(argv[2]);
    int n = atoi(argv[3]);
	int repeats = 1;
	if(argc == 5){
		repeats = atoi(argv[4]);
	}
    
    float *A, *B, *C;
    
    A = new float[m * p];
    B = new float[p * n];
    C = new float[m * n];
        
	double avgtime = 0.0;
	for(int r = 0; r < repeats; r++){ 	
	    initMatrixWithRV(A, m, p);     
	    initMatrixWithRV(B, p, n);
	    initMatrixWithZero(C, m, n);
        
	    clock_t start_time = clock();    

	    matrix_multiply(A, B, C, m, p, n);
    
	    // end timing
	    clock_t end_time = clock();	
		avgtime += (end_time - start_time) / CLOCKS_PER_SEC;
	}
	avgtime = avgtime / repeats;

	std::cout << "m = " << m << ", time cost = " << avgtime << " seconds." << std::endl;
     
    delete[] A;
    delete[] B;
    delete[] C;
    	
    return 0;
}

