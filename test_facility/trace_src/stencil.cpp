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
    #define DIM_SIZE 64
#endif
#ifdef SMALL_DATASET
    #define DIM_SIZE 1024
#endif
#ifdef MEDIUM_DATASET
    #define DIM_SIZE 2048
#endif
#ifdef EXTRALARGE_DATASET
    #define DIM_SIZE 4096
#endif
#ifdef EXTRALARGE_DATASET
    #define DIM_SIZE 8192
#endif


bool varify(double * b, double * a) {
    for (int i = 0; i < DIM_SIZE; i++) {
        for (int j = 0; j < DIM_SIZE; j++) {
            if (b[(i+1) * (DIM_SIZE) + j+1] != a[(i+1)*(DIM_SIZE)+j+1] + a[(i+1)*(DIM_SIZE)+j + 2] + a[(i+1)*(DIM_SIZE)+j] + a[(i)*(DIM_SIZE)+j + 1] + a[(i+2)*(DIM_SIZE)+j+1]) {
                return false;
            }
        }
    }
    return true;
}

void stencil_trace(double *a, double *b, unsigned int dim_size) {

    for (int i = 0; i < DIM_SIZE; i++) {
        for (int j = 0; j < DIM_SIZE; j++) {
            b[(DIM_SIZE)*(i+1)+j+1] = a[(DIM_SIZE)*(i+1)+j+1] 
                        + a[(DIM_SIZE)*(i+1)+j+2] 
                        + a[(DIM_SIZE)*(i+1)+j]
                        + a[(DIM_SIZE)*i+j+1]
                        + a[(DIM_SIZE)*(i+2)+j+1];
            rtTmpAccess(i * (DIM_SIZE) + j + 1, {i, j}, "A[i][j]");
            rtTmpAccess(i * (DIM_SIZE) + j + 2, {i, j}, "A[i][j+1]");
            rtTmpAccess(i * (DIM_SIZE) + j, {i, j}, "A[i][j-1]");
            rtTmpAccess( (i) * (DIM_SIZE) + j + 1, {i, j}, "A[i-1][j]");
            rtTmpAccess( (i+2) * (DIM_SIZE) + j + 1, {i, j}, "A[i+1][j]");
            rtTmpAccess(i * (DIM_SIZE) + j + 1 + (DIM_SIZE) * (DIM_SIZE), {i, j}, "B[i][j]");
/*
            rtTmpAccess((i+1) * (DIM_SIZE) + (j+1));
            rtTmpAccess((i+1) * (DIM_SIZE) + (j+2));
            rtTmpAccess((i+1) * (DIM_SIZE) + (j));
            rtTmpAccess( i * (DIM_SIZE) + (j+1));
            rtTmpAccess( (i+2) * (DIM_SIZE) + (j+1));
            rtTmpAccess(i * (DIM_SIZE) + (j+1) + (DIM_SIZE)* (DIM_SIZE));
            */
        }
    }
    return;
}



int main() {
    double* a = (double*)malloc( (DIM_SIZE+2)* (DIM_SIZE+2)*sizeof(double));
    double* b = (double*)malloc( (DIM_SIZE+2)* (DIM_SIZE+2)*sizeof(double));

    for (int i = 0; i < (DIM_SIZE+2) * (DIM_SIZE+2); i++) {
            a[i] = i % 256;
    }

#ifdef RD
    InitRD();
#endif

#ifdef PAPI_TIMER
    // Get starting timepoint
    PAPI_timer_init();
    PAPI_timer_start();
#endif
    stencil_trace(a, b, DIM_SIZE);
    
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

    if (varify(b, a)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }

    free(a);
    free(b);

    return 0;
}

