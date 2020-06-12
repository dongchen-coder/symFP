#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif


#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
    #define STANDARD_DATASET
#endif
#ifdef MINI_DATASET
    #define DIM_SIZE 64
#endif
#ifdef SMALL_DATASET
    #define DIM_SIZE 1024
#endif
#ifdef STANDARD_DATASET
    #define DIM_SIZE 2048
#endif
#ifdef EXTRALARGE_DATASET
    #define DIM_SIZE 4096
#endif
#ifdef EXTRALARGE_DATASET
    #define DIM_SIZE 8192
#endif


bool varify(double * b, double * a) {
    for (int i = 1; i < DIM_SIZE+1; i++) {
        for (int j = 1; j < DIM_SIZE+1; j++) {
            if (b[i * (DIM_SIZE+2) + j] != a[i*(DIM_SIZE+2)+j] + a[i*(DIM_SIZE+2)+j + 1] + a[i*(DIM_SIZE+2)+j - 1] + a[(i-1)*(DIM_SIZE+2)+j] + a[(i+1)*(DIM_SIZE+2)+j]) {
                return false;
            }
        }
    }
    return true;
}

void stencil_trace(double *a, double *b, unsigned int dim_size) {
    for (int i = 1; i < dim_size+1; i++) {
        for (int j = 1; j < dim_size+1; j++) {
            b[i* (DIM_SIZE + 2) +j] =  a[i* (DIM_SIZE + 2)+j] + a[i* (DIM_SIZE + 2)+j + 1] + a[i* (DIM_SIZE + 2)+j - 1] + a[(i-1)* (DIM_SIZE + 2) +j] + a[(i+1)* (DIM_SIZE + 2) +j];

            rtTmpAccess(i * (DIM_SIZE + 2) + j, {i, j}, "A[i][j]");
            rtTmpAccess(i * (DIM_SIZE + 2) + j + 1, {i, j}, "A[i][j+1]");
            rtTmpAccess(i * (DIM_SIZE + 2) + j - 1, {i, j}, "A[i][j-1]");
            rtTmpAccess( (i-1) * (DIM_SIZE + 2) + j, {i, j}, "A[i-1][j]");
            rtTmpAccess( (i+1) * (DIM_SIZE + 2) + j, {i, j}, "A[i+1][j]");
            rtTmpAccess(i * (DIM_SIZE + 2) + j + (DIM_SIZE + 2)* (DIM_SIZE + 2), {i, j}, "B[i][j]");
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
    
    stencil_trace(a, b, DIM_SIZE);
    
#ifdef PROFILE_RT
    dumpRtTmp();
    RTtoMR_AET();
    dumpMR();
    dumpStat();
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

