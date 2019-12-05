#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

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


#define A_OFFSET 0
#define X1_OFFSET N * N 
#define X2_OFFSET N * N + N
#define Y1_OFFSET N * N + N + N
#define Y2_OFFSET N * N + N + N + N

void runMvt_trace( double* a, double* x1, double* x2, double* y1, double* y2)
{
    int i, j;
    
    for (i=0; i< N; i++)
    {
        for (j=0; j < N; j++)
        {
            //x1[i] = x1[i] + a[i][j] * y1[j];
            x1[i] = x1[i] + a[i * N + j] * y1[j];
		
			rtTmpAccess(X1_OFFSET + i);
			rtTmpAccess(A_OFFSET + i * N + j);
			rtTmpAccess(Y1_OFFSET + j);
			rtTmpAccess(X1_OFFSET + i);			

        }
    }
    
    for (i=0; i < N; i++)
    {
        for (j=0; j< N; j++)
        {
            //x2[i] = x2[i] + a[j][i] * y2[j];
            x2[i] = x2[i] + a[j * N + i] * y2[j];

			rtTmpAccess(X2_OFFSET + i);
			rtTmpAccess(A_OFFSET + j * N + i);
			rtTmpAccess(Y2_OFFSET + j);
			rtTmpAccess(X2_OFFSET + i);		

        }
    }
    return;
}

int main(int argc, char const *argv[])
{
    double* a = (double*)malloc( (N*N) * sizeof(double));
    double* y1 = (double*)malloc( N * sizeof(double));
    double* x1 = (double*)malloc( N * sizeof(double));
    double* y2 = (double*)malloc( N * sizeof(double));
    double* x2 = (double*)malloc( N * sizeof(double));
    

    for (int i = 0; i < N; ++i)
    {
        y1[i] = i % 256;
        y2[i] = i / 100;
    }

    for (int i = 0; i < N*N; ++i)
    {
        a[i] = i / 10;
    }

#ifdef RD
    InitRD();
#endif
    
    runMvt_trace( a, x1, x2, y1, y2);
    
#ifdef PROFILE_RT
    dumpRtTmp();
    RTtoMR_AET();
    dumpMR();
#endif

#ifdef RD
    FiniRD();
#endif

    free(a);
    free(y1);
    free(x1);
    free(y2);
    free(x2);
    
    return 0;
}
