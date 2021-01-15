//void atax_cpu(int nx, int ny, DATA_TYPE POLYBENCH_2D(A,NX,NY,nx,ny), DATA_TYPE POLYBENCH_1D(x,NY,ny), DATA_TYPE POLYBENCH_1D(y,NY,ny), DATA_TYPE POLYBENCH_1D(tmp,NX,nx))

#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif
#ifdef MINI_DATASET
    #define M 32
    #define N 32
#endif
#ifdef SMALL_DATASET
    #define M 1024
    #define N 1024
#endif
#ifdef MEDIUM_DATASET
    #define M 2048
    #define N 2048
#endif
#ifdef LARGE_DATASET
    #define M 8192
    #define N 8192
#endif
#ifdef EXTRALARGE_DATASET
    #define M 100000
    define N 100000
#endif

#else
    #define M 4
    #define N 4
#endif
void atax(int nx, int ny, double* A, double* x, double* y, double* tmp)
{
    int i,j;

    for (i= 0; i < N; i++)
    {
        y[i] = 0;
    }
    
    for (i = 0; i < M; i++)
    {   
        tmp[i] = 0;

  		for (j = 0; j < N; j++)
        {
            tmp[i] = tmp[i] + A[i * N + j] * x[j];
        }
    }
    
    for (i = 0; i < N; i++) 
    {
  		for (j = 0; j < M; j++)
        {
            y[i] = y[i] + A[j * N + i] * tmp[j];
        }
    }
    return;
}
