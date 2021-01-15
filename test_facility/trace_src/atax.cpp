#include "data_size.h"

#ifdef PROFILE_RT
    #include "rt.h"
#endif

#ifdef RD
    #include "reda-spatial.h"
#endif

#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET) && !defined(MEDIUM_DATASET)
    #define LARGE_DATASET
#endif
#ifdef MINI_DATASET
    #define NX 32
    #define NY 32
#endif
#ifdef SMALL_DATASET
    #define NX 1024
    #define NY 1024
#endif
#ifdef MEDIUM_DATASET
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

#define A_OFFSET 0
#define X_OFFSET NX * NY
#define Y_OFFSET NX * NY + NY
#define TMP_OFFSET NX * NY + NY + NY
void atax_cpu_trace(double* A, double* x, double* y, double* tmp, unsigned int dim_size_x, unsigned int dim_size_y) {
    
    int i,j;

    for (i= 0; i < NY; i++)
    {
        y[i] = 0;
        rtTmpAccess(Y_OFFSET + i);
    }

    for (i = 0; i < NX; i++)
    {
            tmp[i] = 0;
            rtTmpAccess(TMP_OFFSET + i);

            for (j = 0; j < NY; j++)
            {
                tmp[i] = tmp[i] + A[i * NY + j] * x[j];
                rtTmpAccess(TMP_OFFSET + i);
                rtTmpAccess(A_OFFSET + i * NY + j);
                rtTmpAccess(X_OFFSET + j);
                rtTmpAccess(TMP_OFFSET + i);
            }

            for (j = 0; j < NY; j++)
            {
                y[j] = y[j] + A[i * NY + j] * tmp[i];
                rtTmpAccess(Y_OFFSET + j);
                rtTmpAccess(A_OFFSET + i * NY + j);
                rtTmpAccess(TMP_OFFSET + i);
                rtTmpAccess(Y_OFFSET + j);
            }
    }


    return;
}


int main() {

    double* A = (double*)malloc( (NX * NY) * sizeof(double) );
    double* x = (double*)malloc( NX * sizeof(double) );
    double* y = (double*)malloc( NY * sizeof(double) );
    double* tmp = (double*)malloc( NX * sizeof(double) );

    for (int i = 0; i < (NX * NY); i++) {
        A[i] = i % 256;
    }

    for (int i = 0; i < NX; i++) {
        x[i] = i % 256;
    }

#ifdef RD
    InitRD();
#endif

#ifdef PAPI_TIMER
    // Get starting timepoint
    PAPI_timer_init();
    PAPI_timer_start();
#endif    
    atax_cpu_trace(A, x, y, tmp, NX, NY);
    
#ifdef PROFILE_RT
    RTtoMR_AET();
#endif 
#ifdef PAPI_TIMER
    // Get ending timepoint
    PAPI_timer_end();
    PAPI_timer_print();
#endif
#ifdef PROFILE_RT
    dumpRtTmp();
    dumpMR();
#endif
    
#ifdef RD
    FiniRD();
#endif

    free(A);
    free(x);
    free(y);
    free(tmp);

    return 0;
}

