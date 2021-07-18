//void covariance(int m, int n, DATA_TYPE POLYBENCH_2D(data,M,N,m,n), DATA_TYPE POLYBENCH_2D(symmat,M,M,m,m), DATA_TYPE POLYBENCH_1D(mean,M,m))

/* Default to LARGE_DATASET. */
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif

# if !defined(M) && !defined(N)
/* Define sample dataset sizes. */
#  ifdef MINI_DATASET
#   define M 16
#   define N 16
#  endif

#  ifdef SMALL_DATASET
#   define M 128
#   define N 128
#  endif

#  ifdef MEDIUM_DATASET
#   define M 2048
#   define N 2048
#  endif

#  ifdef LARGE_DATASET
#   define M 4096
#   define N 4096
#  endif

#  ifdef EXTRALARGE_DATASET
#   define M 8192
#   define N 8192
#  endif


#endif /* !(M N) */


double float_n= 3214212.01;
double eps=  0.005;

void covariance(int m, int n, double* data, double* cov, double* mean)
{
    int i, j, k;
    
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
    for (i = 0; i < M; i++)
    {
        for (j = i; j < M; j++)
        {
            cov[i * M + j] = 0.0;
            for (k = 0; k < N; k++)
            {
                cov[i * M + j] += data[k * M + i] * data[k * M + j];
            }
            cov[j * M + i] = cov[i * M + j];
        }
    }
}
