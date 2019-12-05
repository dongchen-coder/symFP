#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
    #define STANDARD_DATASET
#endif

/* Define the possible dataset sizes. */
#  ifdef MINI_DATASET
#   define NI 32
#   define NJ 32
#  endif

#  ifdef SMALL_DATASET
#   define NI 128
#   define NJ 128
#  endif

#  ifdef STANDARD_DATASET /* Default if unspecified. */
#   define NI 1024
#   define NJ 1024
#  endif

#  ifdef LARGE_DATASET
#   define NI 2048
#   define NJ 2048
#  endif

#  ifdef EXTRALARGE_DATASET
#   define NI 4096
#   define NJ 4096
#  endif

#define A_OFFSET 0
#define C_OFFSET NI * NJ

void syrk_trace(double alpha, double beta, double* A, double* C)
{
    int i, j, k;
    
    /*  C := alpha*A*A' + beta*C */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            C[i * NI + j] = C[i * NI + j] * beta;
			rtTmpAccess(C_OFFSET + i * NI + j);
			rtTmpAccess(C_OFFSET + i * NI + j);
        }
    }
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            for (k = 0; k < NJ; k++)
            {
                C[i * NI + j] = C[i * NI + j] + alpha * A[i * NJ + k] * A[j * NJ + k];
				rtTmpAccess(A_OFFSET + i * NJ + k);
				rtTmpAccess(A_OFFSET + j * NJ + k);
				rtTmpAccess(C_OFFSET + i * NI + j);
				rtTmpAccess(C_OFFSET + i * NI + j);
            }
        }
    }
    return;
}

int main(int argc, char const *argv[])
{
    double* A = (double*)malloc( (NI*NJ) * sizeof(double));
    double* C = (double*)malloc( (NI*NI) * sizeof(double));

    for (int i = 0; i < NI * NJ; ++i)
    {
        A[i] = i / 10;
        C[i] = 1.0;
    }

    double alpha = 0.0;
    double beta = 1.5;

#ifdef RD
    InitRD();
#endif
    
    syrk_trace(alpha, beta, A, C);
    
#ifdef PROFILE_RT
    dumpRtTmp();
    RTtoMR_AET();
    dumpMR();
#endif
    
#ifdef RD
    FiniRD();
#endif

    free(A);
    free(C);

    return 0;
}
