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
#define S_OFFSET NX * NY
#define Q_OFFSET NX * NY + NY
#define R_OFFSET NX * NY + NY + NX
#define P_OFFSET NX * NY + NY + NX + NX


#define A_OFFSET 0
#define S_OFFSET NX * NY
#define Q_OFFSET NX * NY + NY
#define R_OFFSET NX * NY + NY + NX
#define P_OFFSET NX * NY + NY + NX + NX


void bicg_cpu_trace(double* A, double* r, double* s, double* p, double* q, unsigned int nx, int ny) {
 
    int i,j;

    for (i = 0; i < NY; i++)
    {
        s[i] = 0.0;
        rtTmpAccess(S_OFFSET + i);
    }

    for (i = 0; i < NX; i++)
    {
        q[i] = 0.0;
        rtTmpAccess(Q_OFFSET + i);
        for (j = 0; j < NY; j++)
        {
            s[j] = s[j] + r[i] * A[i * NY + j];
            q[i] = q[i] + A[i * NY + j] * p[j];
            rtTmpAccess(S_OFFSET + j);
            rtTmpAccess(R_OFFSET + i);
            rtTmpAccess(A_OFFSET + i * NY + j);
            rtTmpAccess(S_OFFSET + j);
            rtTmpAccess(Q_OFFSET + i);
            rtTmpAccess(A_OFFSET + i * NY + j);
            rtTmpAccess(P_OFFSET + j);
            rtTmpAccess(Q_OFFSET + i);
        }
    }

    return;
}
 


int main() {
    double* A = (double*)malloc( (NX*NY)*sizeof(double));
    double* r = (double*)malloc( NX*sizeof(double));
    double* s = (double*)malloc( NY*sizeof(double));
    double* q = (double*)malloc( NX*sizeof(double));
    double* p = (double*)malloc( NY*sizeof(double));

    for (int i = 0; i < NX; i++) {
        r[i] = i % 256;
        q[i] = 0.0;
    }

    for (int i = 0; i< NY; i++) {
        p[i] = i % 256;
        s[i] = 0.0;
    }

    for (int i = 0; i< NY*NX; i++) {
        A[i] = i % 128;
    }
    
#ifdef RD
    InitRD();
#endif

#ifdef PAPI_TIMER
    // Get starting timepoint
    PAPI_timer_init();
    PAPI_timer_start();
#endif
    
    bicg_cpu_trace(A, r, s, p ,q, NX, NY);
    
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
    free(r);
    free(s);
    free(q);
    free(p);

    return 0;
}

