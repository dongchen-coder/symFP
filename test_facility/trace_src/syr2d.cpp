#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

#ifdef ORG
	#define N 1024
	#define M 1024
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define A_OFFSET 0
#define B_OFFSET M * M
#define C_OFFSET M * M + N * M

void syr2k_trace(double* A, double* B, double* C, double alpha, double beta) {

	int i, j, k;

	for (i = 0; i < N; i++) {
		for (j = 0; j <= i; j++) {
			C[i * N + j] *= beta;
			rtTmpAccess(C_OFFSET + i * N + j);
			rtTmpAccess(C_OFFSET + i * N + j);
		}
		for (k = 0; k < M; k++) {
			for (j = 0; j <= i; j++) {
				C[i * N + j] += A[j * M + k]*alpha*B[i * M + k] + B[j * M + k]*alpha*A[i * N + k];
				rtTmpAccess(A_OFFSET + j * M + k);
				rtTmpAccess(B_OFFSET + i * M + k);
				rtTmpAccess(B_OFFSET + j * M + k);
				rtTmpAccess(A_OFFSET + i * N + k);
				rtTmpAccess(C_OFFSET + i * N + j);
				rtTmpAccess(C_OFFSET + i * N + j);
			}
		}
	}
}

int main() {

	double* A = (double *)malloc(M * M * sizeof(double));
	double* B = (double *)malloc(N * M * sizeof(double));
	double* C = (double *)malloc(N * N * sizeof(double));
	double alpha = 0.2;
	double beta = 0.8;

#ifdef RD
    InitRD();
#endif
    
	syr2k_trace(A, B, C, alpha, beta);

#ifdef PROFILE_RT
    dumpRtTmp();
    RTtoMR_AET();
    dumpMR();
#endif
    
#ifdef RD
    FiniRD();
#endif
    
	return 0;
}

