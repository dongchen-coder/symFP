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
#define u1_OFFSET N * N  
#define v1_OFFSET N * N + N
#define u2_OFFSET N * N + N + N
#define v2_OFFSET N * N + N + N + N
#define W_OFFSET N * N + N + N + N + N 
#define X_OFFSET N * N + N + N + N + N + N
#define Y_OFFSET N * N + N + N + N + N + N + N
#define Z_OFFSET N * N + N + N + N + N + N + N + N

void gemver_trace(double alpha, double beta, double* A, double* u1, double* v1, double* u2, double* v2, double* w, double* x, double* y, double* z)
{
	int i,j;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i * N + j] = A[i * N + j] + u1[i] * v1[j] + u2[i] * v2[j];
			rtTmpAccess(A_OFFSET + i * N + j);
			rtTmpAccess(u1_OFFSET + i);
			rtTmpAccess(v1_OFFSET + j);
			rtTmpAccess(u2_OFFSET + i);
			rtTmpAccess(v2_OFFSET + j);
			rtTmpAccess(A_OFFSET + i * N + j); 
       }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            x[i] = x[i] + beta * A[j * N + i] * y[j];
			rtTmpAccess(X_OFFSET + i);
			rtTmpAccess(A_OFFSET + j * N + i);
			rtTmpAccess(Y_OFFSET + j);
			rtTmpAccess(X_OFFSET + i);
        }
    }

    for (i = 0; i < N; i++)
    {
        x[i] = x[i] + z[i];
		rtTmpAccess(X_OFFSET + i);
		rtTmpAccess(Z_OFFSET + i);
		rtTmpAccess(X_OFFSET + i);
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            w[i] = w[i] +  alpha * A[i * N + j] * x[j];
			rtTmpAccess(W_OFFSET + i);
			rtTmpAccess(A_OFFSET + i * N + j);
			rtTmpAccess(X_OFFSET + j);
			rtTmpAccess(W_OFFSET + i);
        }
    }
}

int main(int argc, char const *argv[])
{
    double* u1 = (double*)malloc( N * sizeof(double));
    double* v1 = (double*)malloc( N * sizeof(double));
    double* u2 = (double*)malloc( N * sizeof(double));
    double* v2 = (double*)malloc( N * sizeof(double));
    double* A = (double*)malloc( (N*N) * sizeof(double));
    double* y = (double*)malloc( N * sizeof(double));
    double* x = (double*)malloc( N * sizeof(double));
    double* z = (double*)malloc( N * sizeof(double));
    double* w = (double*)malloc( N * sizeof(double));

    for (int i = 0; i < N; ++i)
    {
        u1[i] = i % 256;
        u2[i] = i % 32;
        v1[i] = i % 2;
        v2[i] = i % 512;
        y[i] = i % 4;
        z[i] = i % 64;
        x[i] = 0.0;
        w[i] = 0.0;
    }

    for (int i = 0; i < N*N; ++i)
    {
        A[i] = 0.0;
    }

    double alpha = 1.0;
    double beta = 1.5;

#ifdef RD
    InitRD();
#endif
    
    gemver_trace(alpha, beta, A, u1, v1, u2, v2, w, x, y, z);
    
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
