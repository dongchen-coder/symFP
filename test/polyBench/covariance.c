//void covariance(int m, int n, DATA_TYPE POLYBENCH_2D(data,M,N,m,n), DATA_TYPE POLYBENCH_2D(symmat,M,M,m,m), DATA_TYPE POLYBENCH_1D(mean,M,m))

#define M 1024
#define N 1024

double float_n= 3214212.01;
double eps=  0.005;

void covariance(int m, int n, double* data, double* symmat, double* mean)
{
    int i, j, j1,j2;
    
    /* Determine mean of column vectors of input data matrix */
    for (j = 0; j < M; j++)
    {
        mean[j] = 0.0;
        for (i = 0; i < N; i++)
        {
            mean[j] += data[i * M + j];
        }
        mean[j] /= float_n;
    }
    
    /* Center the column vectors. */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            data[i * M + j] -= mean[j];
        }
    }
    
    /* Calculate the m * m covariance matrix. */
    for (j1 = 0; j1 < M; j1++)
    {
        for (j2 = j1; j2 < M; j2++)
        {
            symmat[j1 * M + j2] = 0.0;
            for (i = 0; i < N; i++)
            {
                symmat[j1 * M + j2] += data[i * M + j1] * data[i * M + j2];
            }
            symmat[j2 * M + j1] = symmat[j1 * M + j2];
        }
    }
}
