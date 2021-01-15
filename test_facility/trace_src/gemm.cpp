#include "data_size.h"

#ifdef PROFILE_RT
    #include "rt.h"
#endif

#ifdef RD
    #include "reda-spatial.h"
#endif

#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET) && !defined(MEDIUM_DATASET)
    #define LARGE_DATASET
#endif
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
#ifdef MEDIUM_DATASET
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

#define A_OFFSET 0
#define B_OFFSET NI * NK
#define C_OFFSET NI * NK + NK * NJ

void gemm_trace(double alpha, double beta, double* A, double* B, double* C) {
    int i,j,k;
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            C[i * NJ + j] *= beta;
            rtTmpAccess(C_OFFSET + i * NJ + j); 
            rtTmpAccess(C_OFFSET + i * NJ + j);           

            for (k = 0; k < NK; k++)
            {
                C[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
                rtTmpAccess(A_OFFSET + i * NK + k);
                rtTmpAccess(B_OFFSET + k * NJ + j);
                rtTmpAccess(C_OFFSET + i * NJ + j);
                rtTmpAccess(C_OFFSET + i * NJ + j);
            }
        }
    }
    return;
}

int main()
{
    double* A = (double*)malloc(NI * NK * sizeof(double));
    double* B = (double*)malloc(NK * NJ * sizeof(double));
    double* C = (double*)malloc(NI * NJ * sizeof(double));

    for (int i = 0; i < NI*NK; ++i)
    {
        A[i] = i % 256;
    }

    for (int i = 0; i < NK * NJ; ++i)
    {
        B[i] = i % 48;
    }

    double alpha = 1.0;
    double beta = 1.5;

#ifdef RD
    InitRD();
#endif

#ifdef PAPI_TIMER
    // Get starting timepoint
    PAPI_timer_init();
    PAPI_timer_start();
#endif
    
    gemm_trace(alpha, beta, A, B, C);
    
#ifdef PROFILE_RT
    RTtoMR_AET();
#endif 
#ifdef PAPI_TIMER
    // Get ending timepoint
    PAPI_timer_end();
    PAPI_timer_print();
#endif
#ifdef PROFILE_RT
    dumpRtTmp();
    dumpMR();
#endif
    
#ifdef RD
    FiniRD();
#endif

    free(A);
    free(B);
    free(C);

    return 0;
}
