//void gemver(int n, DATA_TYPE alpha, DATA_TYPE beta, DATA_TYPE POLYBENCH_2D(A, N, N, n, n), DATA_TYPE POLYBENCH_1D(u1, N, n), DATA_TYPE POLYBENCH_1D(v1, N, n), DATA_TYPE POLYBENCH_1D(u2, N, n), DATA_TYPE POLYBENCH_1D(v2, N, n), DATA_TYPE POLYBENCH_1D(w, N, n), DATA_TYPE POLYBENCH_1D(x, N, n), DATA_TYPE POLYBENCH_1D(y, N, n), DATA_TYPE POLYBENCH_1D(z, N, n))

#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif    
#ifdef MINI_DATASET
    #define N 32
#endif
#ifdef SMALL_DATASET
    #define N 1024
#endif
#ifdef MEDIUM_DATASET
    #define N 2048
#endif
#ifdef LARGE_DATASET
    #define N 8192
#endif
#ifdef EXTRALARGE_DATASET
    #define N 100000
#endif

#else 
    #define N 4
#endif

void gemver(int n, double alpha, double beta, double* A, double* u1, double* v1, double* u2, double* v2, double* w, double* x, double* y, double* z)
{
    int i,j;
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i * N + j] = A[i * N + j] + u1[i] * v1[j] + u2[i] * v2[j];
        }
    }
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            x[i] = x[i] + beta * A[j * N + i] * y[j];
        }
    }
    
    for (i = 0; i < N; i++)
    {
        x[i] = x[i] + z[i];
    }
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            w[i] = w[i] +  alpha * A[i * N + j] * x[j];
        }
    }
}
