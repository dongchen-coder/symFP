#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

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

		rtTmpAccess(TMP_OFFSET + i);
		rtTmpAccess(Y_OFFSET + i);
	
        for (j = 0; j < N; j++)
        {
            tmp[i] = A[i * N + j] * x[j] + tmp[i];
            y[i] = B[i * N + j] * x[j] + y[i];
        	rtTmpAccess(A_OFFSET + i * N + j);
			rtTmpAccess(X_OFFSET + j);
			rtTmpAccess(TMP_OFFSET + i);
			rtTmpAccess(TMP_OFFSET + i);
			rtTmpAccess(B_OFFSET + i * N + j);
            rtTmpAccess(X_OFFSET + j);
            rtTmpAccess(Y_OFFSET + i);
            rtTmpAccess(Y_OFFSET + i);
		}

        y[i] = alpha * tmp[i] + beta * y[i];

		rtTmpAccess(TMP_OFFSET + i);
		rtTmpAccess(Y_OFFSET + i);
		rtTmpAccess(Y_OFFSET + i);
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

#ifdef RD
    InitRD();
#endif
    
    gesummv_trace(alpha, beta, A, B, tmp, x, y);
    
#ifdef PROFILE_RT
    dumpRtTmp();
    RTtoMR_AET();
    dumpMR();
#endif
    
#ifdef RD
    FiniRD();
#endif

    return 0;
}
