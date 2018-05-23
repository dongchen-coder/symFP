#include "../utility/data_size.h"

#ifdef ORG
    #define NI 256
    #define NJ 256
    #define NL 256
    #define NK 256
    #define NM 256
#elif defined (TX)
    #define NI 512
    #define NJ 256
    #define NL 512
    #define NK 256
    #define NM 256
#elif defined (FX)
    #define NI 1024
    #define NJ 256
    #define NL 1024
    #define NK 256
    #define NM 256
#elif defined (EX)
    #define NI 2048
    #define NJ 256
    #define NL 2048
    #define NK 256
    #define NM 256
#endif


void mm3_cpu(int ni, int nj, int nk, int nl, int nm,
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
