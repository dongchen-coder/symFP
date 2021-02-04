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
    #define N 32
#endif
#ifdef SMALL_DATASET
    #define N 1024
#endif
#ifdef MEDIUM_DATASET
    #define N 4096
#endif
#ifdef LARGE_DATASET
    #define N 8192
#endif
#ifdef EXTRALARGE_DATASET
    #define N 100000
#endif


#define A_OFFSET 0
#define B_OFFSET N * N
#define TMP_OFFSET N * N + N * N
#define X_OFFSET N * N + N * N + N
#define Y_OFFSET N * N + N * N + N + N

void gesummv_trace(double alpha, double beta, double* A, double* B, double* tmp, double* x, double* y)
{
    int i, j;    

    for (i = 0; i < N; i++)
    {
        tmp[i] = 0;
        y[i] = 0;

        // rtTmpAccess(TMP_OFFSET + i);
        // rtTmpAccess(Y_OFFSET + i);
        rtTmpAccess(TMP_OFFSET + i, "tmp_addr0", {i});
        rtTmpAccess(Y_OFFSET + i, "y_addr0", {i});
    
        for (j = 0; j < N; j++)
        {
            tmp[i] = A[i * N + j] * x[j] + tmp[i];
            y[i] = B[i * N + j] * x[j] + y[i];
            // rtTmpAccess(A_OFFSET + i * N + j);
            // rtTmpAccess(X_OFFSET + j);
            // rtTmpAccess(TMP_OFFSET + i);
            // rtTmpAccess(TMP_OFFSET + i);
            // rtTmpAccess(B_OFFSET + i * N + j);
            // rtTmpAccess(X_OFFSET + j);
            // rtTmpAccess(Y_OFFSET + i);
            // rtTmpAccess(Y_OFFSET + i);
            rtTmpAccess(A_OFFSET + i * N + j, "A_addr0", {i, j});
            rtTmpAccess(X_OFFSET + j, "X_addr0", {i, j});
            rtTmpAccess(TMP_OFFSET + i, "tmp_addr1", {i, j});
            rtTmpAccess(TMP_OFFSET + i, "tmp_addr2", {i, j});
            rtTmpAccess(B_OFFSET + i * N + j, "B_addr0", {i, j});
            rtTmpAccess(X_OFFSET + j, "x_addr1", {i, j});
            rtTmpAccess(Y_OFFSET + i, "y_addr1", {i, j});
            rtTmpAccess(Y_OFFSET + i, "y_addr2", {i, j});
        }

        y[i] = alpha * tmp[i] + beta * y[i];

        // rtTmpAccess(TMP_OFFSET + i);
        // rtTmpAccess(Y_OFFSET + i);
        // rtTmpAccess(Y_OFFSET + i);
        rtTmpAccess(TMP_OFFSET + i, "tmp_addr3", {i});
        rtTmpAccess(Y_OFFSET + i, "y_addr3", {i});
        rtTmpAccess(Y_OFFSET + i, "y_addr4", {i});
    }

    return;
}

int main(int argc, char const *argv[])
{
    double* A = (double*)malloc( (N*N) * sizeof(double));
    double* x = (double*)malloc( N * sizeof(double));
    double* tmp = (double*)malloc( N * sizeof(double));
    double* B = (double*)malloc( (N*N) * sizeof(double));
    double* y = (double*)malloc( N * sizeof(double));

    for (int i = 0; i < N; ++i)
    {
        x[i] = i % 256;
    }

    for (int i = 0; i < N*N; ++i)
    {
        A[i] = i / 10;
        B[i] = i / 25;
    }

    double alpha = 1.0;
    double beta = 1.5;

#ifdef RD
    InitRD();
#endif
#ifdef PAPI_TIMER
    // Get starting timepoint
    PAPI_timer_init();
    PAPI_timer_start();
#endif
    
    gesummv_trace(alpha, beta, A, B, tmp, x, y);
    
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
    free(tmp);
    free(x);
    free(B);
    free(y);

    return 0;
}
