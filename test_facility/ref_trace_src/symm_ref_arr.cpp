

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
#define C_OFFSET M * M + M * N

void symm_trace(double* A, double* B, double* C, double alpha, double beta) {
	
	int i, j, k;	
	double temp2;

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++ ) {
			temp2 = 0;
			for (k = 0; k < i; k++) {
				C[k * N + j] += alpha*B[i * N + j] * A[i * M + k];
				temp2 += B[k * N + j] * A[i * M + k];
                rtTmpAccess(C_OFFSET + k * N + j, 0, 0);
                rtTmpAccess(B_OFFSET + i * N + j, 1, 1);
                rtTmpAccess(A_OFFSET + i * M + k, 2, 2);
                rtTmpAccess(B_OFFSET + k * N + j, 3, 1);
                rtTmpAccess(A_OFFSET + i * M + k, 4, 2);
			}
			C[i * N + j] = beta * C[i * N + j] + alpha*B[i * N + j] * A[i * M + i] + alpha * temp2;
            rtTmpAccess(C_OFFSET + i * N + j, 5, 0);
            rtTmpAccess(B_OFFSET + i * N + j, 6, 1);
            rtTmpAccess(A_OFFSET + i * M + i, 7, 2);
            rtTmpAccess(C_OFFSET + i * N + j, 8, 0);
		}
	}
	return;
}

int main() {

	double* A = (double *)malloc(M * M * sizeof(double));
	double* B = (double *)malloc(M * N * sizeof(double));
	double* C = (double *)malloc(M * N * sizeof(double));
	double alpha = 0.2;
	double beta = 0.8;

	symm_trace(A, B, C, alpha, beta);

	dumpSetSize();

	return 0;
}
