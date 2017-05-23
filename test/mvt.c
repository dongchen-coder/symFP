//void runMvt(int n, DATA_TYPE POLYBENCH_2D(a, N, N, n, n), DATA_TYPE POLYBENCH_1D(x1, N, n), DATA_TYPE POLYBENCH_1D(x2, N, n), DATA_TYPE POLYBENCH_1D(y1, N, n), DATA_TYPE POLYBENCH_1D(y2, N, n))

#define N 1024

void runMvt(int n, double* a, double* x1, double* x2, double* y1, double* y2)
{
    int i, j;
    
    for (i=0; i< N; i++)
    {
        for (j=0; j < N; j++)
        {
            //x1[i] = x1[i] + a[i][j] * y1[j];
            x1[i] = x1[i] + a[i * N + j] * y1[j];
        }
    }
    
    for (i=0; i < N; i++)
    {
        for (j=0; j< N; j++)
        {
            //x2[i] = x2[i] + a[j][i] * y2[j];
            x2[i] = x2[i] + a[j * N + i] * y2[j];
        }
    }
}
