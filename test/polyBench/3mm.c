/*void mm3_cpu(int ni, int nj, int nk, int nl, int nm,
             DATA_TYPE POLYBENCH_2D(E,NI,NJ,ni,nj),
             DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
             DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj),
             DATA_TYPE POLYBENCH_2D(F,NJ,NL,nj,nl),
             DATA_TYPE POLYBENCH_2D(C,NJ,NM,nj,nm),
             DATA_TYPE POLYBENCH_2D(D,NM,NL,nm,nl),
             DATA_TYPE POLYBENCH_2D(G,NI,NL,ni,nl))
*/

#define NI 1024
#define NJ 1024
#define NL 1024
#define NK 1024
#define NM 1024

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
