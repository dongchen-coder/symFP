//void correlation(int m, int n, DATA_TYPE POLYBENCH_2D(data, M, N, m, n), DATA_TYPE POLYBENCH_1D(mean, M, m), DATA_TYPE POLYBENCH_1D(stddev, M, m), DATA_TYPE POLYBENCH_2D(symmat, M, N, m, n))

#define M 1024
#define N 1024
#define FLOAT_N 3214212.01
#define EPS 0.005

void correlation(int m, int n, double* data, double* mean, double* stddev, double* symmat)
{
    int i, j, j1, j2;
    
    // Determine mean of column vectors of input data matrix
    for (j = 0; j < M; j++)
   	{
        mean[j] = 0.0;
        
        for (i = 0; i < N; i++)
        {
            mean[j] += data[i * M + j];
        }
        
        mean[j] /= (double)FLOAT_N;
   	}
    
    // Determine standard deviations of column vectors of data matrix.
    for (j = 0; j < M; j++)
   	{
        stddev[j] = 0.0;
        
        for (i = 0; i < N; i++)
        {
            stddev[j] += (data[i * M + j] - mean[j]) * (data[i * M + j] - mean[j]);
        }
        
        stddev[j] /= FLOAT_N;
        stddev[j] = sqrt_of_array_cell(stddev, j);
        stddev[j] = stddev[j] <= EPS ? 1.0 : stddev[j];
    }
    
    // Center and reduce the column vectors.
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            data[i * M + j] -= mean[j];
            data[i * M + j] /= (sqrt(FLOAT_N)*stddev[j]) ;
        }
    }
    
    // Calculate the m * m correlation matrix.
    for (j1 = 0; j1 < M-1; j1++)
    {
        symmat[j1 * M + j1] = 1.0;
        
        for (j2 = j1+1; j2 < M; j2++)
        {
            symmat[j1 * M + j2] = 0.0;
            
            for (i = 0; i < N; i++)
            {
                symmat[j1 * M + j2] += (data[i * M + j1] * data[i * M + j2]);
            }
            
            symmat[j2 * M + j1] = symmat[j1 * M + j2];
        }
    }
    
    symmat[(M-1) * M + M-1] = 1.0;
    
    return;
}

