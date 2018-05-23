
#define N 1024
#define M 1024

void syr2k(double* A, double* B, double* C, double alpha, double beta) {

	int i, j, k;

	for (i = 0; i < N; i++) {
		for (j = 0; j <= i; j++) {
			C[i * N + j] *= beta;
		}
		for (k = 0; k < M; k++) {
			for (j = 0; j <= i; j++) {
				C[i * N + j] += A[j * M + k]*alpha*B[i * M + k] + B[j * M + k]*alpha*A[i * N + k];
			}
		}
	}
}

