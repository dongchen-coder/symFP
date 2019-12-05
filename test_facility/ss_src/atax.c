//void atax_cpu(int nx, int ny, DATA_TYPE POLYBENCH_2D(A,NX,NY,nx,ny), DATA_TYPE POLYBENCH_1D(x,NY,ny), DATA_TYPE POLYBENCH_1D(y,NY,ny), DATA_TYPE POLYBENCH_1D(tmp,NX,nx))

#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define STANDARD_DATASET
# endif
#ifdef MINI_DATASET
    #define NX 32
    #define NY 32
#endif
#ifdef SMALL_DATASET
    #define NX 1024
    #define NY 1024
#endif
#ifdef STANDARD_DATASET
    #define NX 4096
    #define NY 4096
#endif
#ifdef LARGE_DATASET
    #define NX 8192
    #define NY 8192
#endif
#ifdef EXTRALARGE_DATASET
    #define NX 100000
    #define NY 100000
#endif

#else
    #define NX 4
    #define NY 4
#endif
void atax(int nx, int ny, double* A, double* x, double* y, double* tmp)
{
    int i,j;

    // for (i= 0; i < NY; i++)
    // {
    //     y[i] = 0;
    // }
    
    for (i = 0; i < NX; i++)
    {   
        // tmp[i] = 0;

  		for (j = 0; j < NY; j++)
        {
            tmp[i] = tmp[i] + A[i * NY + j] * x[j];
        }
    
  		for (j = 0; j < NY; j++)
        {
            y[j] = y[j] + A[i * NY + j] * tmp[i];
        }
    }
    return;
}
