//void gesummv(int n, DATA_TYPE alpha, DATA_TYPE beta, DATA_TYPE POLYBENCH_2D(A,N,N,n,n), DATA_TYPE POLYBENCH_2D(B,N,N,n,n), DATA_TYPE POLYBENCH_1D(tmp,N,n), DATA_TYPE POLYBENCH_1D(x,N,n), DATA_TYPE POLYBENCH_1D(y,N,n))
#define N 1024

void gesummv(int n, double alpha, double beta, double* A, double* B, double* tmp, double* x, double* y)
{
    int i, j;
    
    for (i = 0; i < N; i++)
    {
        tmp[i] = 0;
        y[i] = 0;
        for (j = 0; j < N; j++)
        {
            tmp[i] = A[i * N + j] * x[j] + tmp[i];
            y[i] = B[i * N + j] * x[j] + y[i];
        }
        
        y[i] = alpha * tmp[i] + beta * y[i];
    }
}
