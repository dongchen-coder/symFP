#include "../utility/rt.h"

#define NI 256
#define NJ 256

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
				rtTmpAccess(C_OFFSET + i * NI + j);
				rtTmpAccess(A_OFFSET + i * NJ + k);
				rtTmpAccess(A_OFFSET + j * NJ + k);
				rtTmpAccess(C_OFFSET + i * NI + j);
            }
        }
    }
    return;
}

int main(int argc, char const *argv[])
{
    double* A = (double*)malloc( (NI*NJ) * sizeof(double));
    double* C = (double*)malloc( (NI*NJ) * sizeof(double));

    for (int i = 0; i < NI * NJ; ++i)
    {
        A[i] = i / 10;
        C[i] = 1.0;
    }

    double alpha = 0.0;
    double beta = 1.5;

    syrk_trace(alpha, beta, A, C);
    dumpRtTmp();
	RTtoMR_AET();
    dumpMR();

    return 0;
}
