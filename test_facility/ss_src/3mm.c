#ifndef DEBUG
#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
    #define STANDARD_DATASET
#endif
#ifdef MINI_DATASET
    #define NI 32
    #define NJ 32
    #define NK 32
    #define NL 32
    #define NM 32
#endif
#ifdef SMALL_DATASET
    #define NI 128
    #define NJ 128
    #define NL 128
    #define NK 128
    #define NM 128
#endif 
#ifdef STANDARD_DATASET
    #define NI 1024
    #define NJ 1024
    #define NK 1024
    #define NL 1024
    #define NM 1024
#endif
#ifdef LARGE_DATASET
    #define NI 2048
    #define NJ 2048
    #define NL 2048
    #define NK 2048
    #define NM 2048
#endif
#ifdef EXTRALARGE_DATASET
    #define NI 4096
    #define NJ 4096
    #define NL 4096
    #define NK 4096
    #define NM 4096
#endif

#else
    #define NI 4
    #define NJ 4
    #define NL 4
    #define NK 4
    #define NM 4
#endif
void mm2(int ni, int nj, int nk, int nl, int nm,
             double * E, double* A, double* B, double* F, double* C, double* D, double* G)
{
    int i, j, k;
    
    /* E := A*B */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            E[i * NJ + j] = 0;
            for (k = 0; k < NK; ++k)
            {
                E[i * NJ + j] += A[i * NK + k] * B[k * NJ + j];
            }
        }
    }
    
    /* F := C*D */
    for (i = 0; i < NJ; i++)
    {
        for (j = 0; j < NL; j++)
        {
            F[i * NL + j] = 0;
            for (k = 0; k < NM; ++k)
            {
                F[i * NL + j] += C[i * NM + k] * D[k * NL + j];
            }
        }
    }
    
    /* G := E*F */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NL; j++)
        {
            G[i * NL + j] = 0;
            for (k = 0; k < NJ; ++k)
            {
                G[i * NL + j] += E[i * NJ + k] * F[k * NL + j];
            }
        }
    }
}
