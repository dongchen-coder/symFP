//void runMvt(int n, DATA_TYPE POLYBENCH_2D(a, N, N, n, n), DATA_TYPE POLYBENCH_1D(x1, N, n), DATA_TYPE POLYBENCH_1D(x2, N, n), DATA_TYPE POLYBENCH_1D(y1, N, n), DATA_TYPE POLYBENCH_1D(y2, N, n))

#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define STANDARD_DATASET
# endif    
#ifdef MINI_DATASET
    #define N 32
#endif
#ifdef SMALL_DATASET
    #define N 1024
#endif
#ifdef STANDARD_DATASET
    #define N 4096
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

void mvt(int n, double* a, double* x1, double* x2, double* y1, double* y2)
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


