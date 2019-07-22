
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define NI 256
	#define NJ 256
	#define NK 256
	#define NL 256
#elif defined (TX)
#elif defined (FX)
#elif defined (EX)
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
            rtTmpAccess(TMP_OFFSET + i * NJ + j, 0, 0);
            for (k = 0; k < NK; ++k) {
                tmp[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
                rtTmpAccess(A_OFFSET + i * NK + k, 1, 1);
                rtTmpAccess(B_OFFSET + k * NJ + j, 2, 2);
                rtTmpAccess(TMP_OFFSET + i * NJ + j, 3, 0);
                rtTmpAccess(TMP_OFFSET + i * NJ + j, 4, 0);
            }
        }
    }
    for (i = 0; i < NI; i++) {
        for (j = 0; j < NL; j++) {
            D[i * NL + j] *= beta;
            rtTmpAccess(D_OFFSET + i * NL + j, 5, 3);
            rtTmpAccess(D_OFFSET + i * NL + j, 6, 3);
            for (k = 0; k < NJ; ++k) {
                D[i * NL + j] += tmp[i * NJ + k] * C[k * NL + j];
                rtTmpAccess(TMP_OFFSET + i * NJ + k, 7, 0);
                rtTmpAccess(C_OFFSET + k * NL + j, 8, 4);
                rtTmpAccess(D_OFFSET + i * NL + j, 9, 3);
                rtTmpAccess(D_OFFSET + i * NL + j, 10, 3);
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
	
	mm2_trace(tmp, A, B, C, D, alpha, beta);

	OSL_ref(0);

	return 0;
}

