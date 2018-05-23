//void gemm(int ni, int nj, int nk, DATA_TYPE alpha, DATA_TYPE beta, DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk), DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj), DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj))

#include "../utility/data_size.h"


#ifdef ORG
    #define NI 256
    #define NJ 256
    #define NK 256
#elif defined(TX)
    #define NI 256
    #define NJ 256
    #define NK 512
#elif defined(FX)
    #define NI 256
    #define NJ 512
    #define NK 512
#elif defined(EX)
    #define NI 512
    #define NJ 512
    #define NK 512
#endif


void gemm(int ni, int nj, int nk, double alpha, double beta, double* A, double* B, double* C)
{
    int i,j,k;
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            C[i * NJ + j] *= beta;
            
            for (k = 0; k < NK; k++)
            {
                C[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
            }
        }
    }
}
