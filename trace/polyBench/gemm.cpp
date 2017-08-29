#include "../utility/rt.h"

#define NI 256
#define NJ 256
#define NK 256

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
				rtTmpAccess(C_OFFSET + i * NJ + j);
				rtTmpAccess(A_OFFSET + i * NK + k);
				rtTmpAccess(B_OFFSET + k * NJ + j);
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

    gemm_trace(alpha, beta, A, B, C);
    dumpRtTmp();

    return 0;
}
