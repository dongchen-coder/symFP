//void gemm(int ni, int nj, int nk, DATA_TYPE alpha, DATA_TYPE beta, DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk), DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj), DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj))

#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define STANDARD_DATASET
# endif
#ifdef MINI_DATASET
    #define NI 32
    #define NJ 32
    #define NK 32
#endif
#ifdef SMALL_DATASET
    #define NI 256
    #define NJ 256
    #define NK 256
#endif
#ifdef STANDARD_DATASET
    #define NI 1024
    #define NJ 1024
    #define NK 1024
#endif
#ifdef LARGE_DATASET
    #define NI 2048
    #define NJ 2048
    #define NK 2048
#endif
#ifdef EXTRALARGE_DATASET
    #define NI 4096
    #define NJ 4096
    #define NK 4096
#endif

#else 
    #define NI 4
    #define NJ 4
    #define NK 4
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
