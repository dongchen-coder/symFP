#include "../utility/rt.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define M 1024
	#define N 1024
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define A_OFFSET 0
#define B_OFFSET M * M

void trmm_trace(double* A, double* B, double alpha) {

	int i, j, k;

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			for (k = i+1; k < M; k++) {
				B[i * N + j] += A[k * M + i] * B[k * N + j];
				rtTmpAccess(A_OFFSET + k * M + i);
				rtTmpAccess(B_OFFSET + k * N + j);
				rtTmpAccess(B_OFFSET + i * N + j);
				rtTmpAccess(B_OFFSET + i * N + j);
			}
			B[i * N + j] = alpha * B[i * N + j];
     		rtTmpAccess(B_OFFSET + i * N + j);
			rtTmpAccess(B_OFFSET + i * N + j);
		}
	}
}

int main() {

	double* A = (double *)malloc(M * M * sizeof(double));
	double* B = (double *)malloc(M * N * sizeof(double));
	double alpha = 0.2;

	trmm_trace(A, B, alpha);

	dumpRtTmp();
    RTtoMR_AET();
    dumpMR();

	return 0;
}
