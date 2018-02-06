
#define M 1024
#define N 1024

void trmm(double* A, double* B, double alpha) {

	int i, j, k;

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			for (k = i+1; k < M; k++) {
				B[i * N + j] += A[k * M + i] * B[k *N + j];
			}
			B[i * N + j] = alpha * B[i * N + j];
     	}
	}
}
