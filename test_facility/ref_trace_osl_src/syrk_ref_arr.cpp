
#include "../utility/reference_lease.h"
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
            rtTmpAccess(C_OFFSET + i * NI + j, 0, 0);
            rtTmpAccess(C_OFFSET + i * NI + j, 1, 0);
        }
    }
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            for (k = 0; k < NJ; k++)
            {
                C[i * NI + j] = C[i * NI + j] + alpha * A[i * NJ + k] * A[j * NJ + k];
                rtTmpAccess(A_OFFSET + i * NJ + k, 2, 1);
                rtTmpAccess(A_OFFSET + j * NJ + k, 3, 1);
                rtTmpAccess(C_OFFSET + i * NI + j, 4, 0);
                rtTmpAccess(C_OFFSET + i * NI + j, 5, 0);
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

    syrk_trace(alpha, beta, A, C);
    
    OSL_ref(0);

    return 0;
}
