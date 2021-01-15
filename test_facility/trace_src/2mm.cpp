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
    #define NI 32
    #define NJ 32
    #define NK 32
    #define NL 32
#endif
#ifdef SMALL_DATASET
    #define NI 128
    #define NJ 128
    #define NL 128
    #define NK 128
#endif
#ifdef MEDIUM_DATASET
    #define NI 1024
    #define NJ 1024
    #define NK 1024
    #define NL 1024
#endif
#ifdef LARGE_DATASET
    #define NI 2048
    #define NJ 2048
    #define NL 2048
    #define NK 2048
#endif
#ifdef EXTRALARGE_DATASET
    #define NI 4096
    #define NJ 4096
    #define NL 4096
    #define NK 4096
#endif


#define TMP_OFFSET 0
#define A_OFFSET NI * NJ
#define B_OFFSET NI * NJ + NI * NK
#define D_OFFSET NI * NJ + NI * NK + NJ * NK
#define C_OFFSET NI * NJ + NI * NK + NJ * NK + NI * NL



void mm2_trace(double* tmp, double* A, double* B, double* C, double* D, double alpha, double beta) {

    int i, j, k;

    for (i = 0; i < NI; i++) {
        for (j = 0; j < NJ; j++) {
            tmp[i * NJ + j] = 0.0;
            rtTmpAccess(TMP_OFFSET + i * NJ + j);
            for (k = 0; k < NK; ++k) {
                tmp[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
                rtTmpAccess(A_OFFSET + i * NK + k);
                rtTmpAccess(B_OFFSET + k * NJ + j);
                rtTmpAccess(TMP_OFFSET + i * NJ + j);
                rtTmpAccess(TMP_OFFSET + i * NJ + j);
            }
        }
    }
    for (i = 0; i < NI; i++) {
        for (j = 0; j < NL; j++) {
            D[i * NL + j] *= beta;
            rtTmpAccess(D_OFFSET + i * NL + j);
            rtTmpAccess(D_OFFSET + i * NL + j);
            for (k = 0; k < NJ; ++k) {
                D[i * NL + j] += tmp[i * NJ + k] * C[k * NL + j];
                rtTmpAccess(TMP_OFFSET + i * NJ + k);
                rtTmpAccess(C_OFFSET + k * NL + j);
                rtTmpAccess(D_OFFSET + i * NL + j);
                rtTmpAccess(D_OFFSET + i * NL + j);
            }
        }
    }
}


int main() {
	
	double* tmp = (double*)malloc( NI * NJ * sizeof(double));
	double* A = (double*)malloc( NI * NK * sizeof(double));
	double* B = (double*)malloc( NK * NJ * sizeof(double));
	double* C = (double*)malloc( NJ * NL * sizeof(double));
	double* D = (double*)malloc( NI * NL * sizeof(double));
	double alpha = 0.1;
	double beta = 0.5;
	
#ifdef RD
	InitRD();
#endif

#ifdef PAPI_TIMER
    // Get starting timepoint
    PAPI_timer_init();
    PAPI_timer_start();
#endif

	mm2_trace(tmp, A, B, C, D, alpha, beta);

#ifdef PROFILE_RT
    RTtoMR_AET();
#endif

#ifdef PROFILE_RT
    dumpRtTmp();
    dumpMR();
#endif


#ifdef RD
    FiniRD();
#endif

    free(A);
    free(B);
    free(C);
    free(D);
    free(tmp);

	return 0;
}


