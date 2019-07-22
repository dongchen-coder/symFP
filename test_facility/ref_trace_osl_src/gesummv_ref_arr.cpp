
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define N 1024
#elif defined(TX)
	#define N 1448
#elif defined(FX)
	#define N 2048
#elif defined(EX)
	#define N 2896
#endif


#define A_OFFSET 0
#define B_OFFSET N * N
#define TMP_OFFSET N * N + N * N
#define X_OFFSET N * N + N * N + N
#define Y_OFFSET N * N + N * N + N + N

void gesummv_trace(double alpha, double beta, double* A, double* B, double* tmp, double* x, double* y)
{
    int i, j;
    
    for (i = 0; i < N; i++)
    {
        tmp[i] = 0;
        y[i] = 0;
        
        rtTmpAccess(TMP_OFFSET + i, 0, 0);
        rtTmpAccess(Y_OFFSET + i, 1, 1);
        
        for (j = 0; j < N; j++)
        {
            tmp[i] = A[i * N + j] * x[j] + tmp[i];
            y[i] = B[i * N + j] * x[j] + y[i];
            rtTmpAccess(A_OFFSET + i * N + j, 2, 2);
            rtTmpAccess(X_OFFSET + j, 3, 3);
            rtTmpAccess(TMP_OFFSET + i, 4, 0);
            rtTmpAccess(TMP_OFFSET + i, 5, 0);
            rtTmpAccess(B_OFFSET + i * N + j, 6, 4);
            rtTmpAccess(X_OFFSET + j, 7, 3);
            rtTmpAccess(Y_OFFSET + i, 8, 1);
            rtTmpAccess(Y_OFFSET + i, 9, 1);
        }
        
        y[i] = alpha * tmp[i] + beta * y[i];
        
        rtTmpAccess(TMP_OFFSET + i, 10, 0);
        rtTmpAccess(Y_OFFSET + i, 11, 1);
        rtTmpAccess(Y_OFFSET + i, 12, 1);
    }
    
    return;
}

int main(int argc, char const *argv[])
{
    double* A = (double*)malloc( (N*N) * sizeof(double));
    double* x = (double*)malloc( N * sizeof(double));
    double* tmp = (double*)malloc( N * sizeof(double));
    double* B = (double*)malloc( (N*N) * sizeof(double));
    double* y = (double*)malloc( N * sizeof(double));

    for (int i = 0; i < N; ++i)
    {
        x[i] = i % 256;
    }

    for (int i = 0; i < N*N; ++i)
    {
        A[i] = i / 10;
        B[i] = i / 25;
    }

    double alpha = 1.0;
    double beta = 1.5;

    gesummv_trace(alpha, beta, A, B, tmp, x, y);
    
    OSL_ref(0);

    return 0;
}
