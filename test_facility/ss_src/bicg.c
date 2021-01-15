//void bicg_cpu(int nx, int ny, DATA_TYPE POLYBENCH_2D(A,NX,NY,nx,ny), DATA_TYPE POLYBENCH_1D(r,NX,nx), DATA_TYPE POLYBENCH_1D(s,NY,ny), DATA_TYPE POLYBENCH_1D(p,NY,ny), DATA_TYPE POLYBENCH_1D(q,NX,nx))

#ifndef DEBUG
#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
    #define LARGE_DATASET
#endif
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
    #define N 100000
#endif

#else 
    #define M 4
    #define N 4
#endif

void bicg_cpu(int nx, int ny, double* A, double* r, double* s, double* p, double* q) {
    int i,j;

    for (i = 0; i < M; i++)
    {
        s[i] = 0.0;
    }

    for (i = 0; i < N; i++)
    {
        q[i] = 0.0;
        for (j = 0; j < M; j++)
        {
            q[i] = q[i] + A[i * M + j] * p[j];
        }
    }

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            s[i] = s[i] + r[j] * A[j * M + i];
        }
    }
}
