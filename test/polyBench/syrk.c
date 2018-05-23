//void syrk(int ni, int nj, DATA_TYPE alpha, DATA_TYPE beta, DATA_TYPE POLYBENCH_2D(A, NI, NJ, ni, nj), DATA_TYPE POLYBENCH_2D(C, NI, NI, ni, ni))

#include "../utility/data_size.h"

#ifdef ORG
    #define NI 256
    #define NJ 256
#elif defined(TX)
    #define NI 362
    #define NJ 256
#elif defined(FX)
    #define NI 512
    #define NJ 256
#elif defined(EX)
    #define NI 724
    #define NJ 256
#endif

void syrk(int ni, int nj, double alpha, double beta, double* A, double* C)
{
    int i, j, k;
    
    /*  C := alpha*A*A' + beta*C */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            C[i * NI + j] *= beta;
        }
    }
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            for (k = 0; k < NJ; k++)
            {
                C[i * NI + j] += alpha * A[i * NJ + k] * A[j * NJ + k];
            }
        }
    }
}
