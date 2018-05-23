
#define N 1024


void cholesky(double* A) {

	int i, j, k;

	for (i = 0; i < N; i++) {
		for (j = 0; j < i; j++) {
			for (k = 0; k < j; k++) {
				A[i * N + j] -= A[i * N + k] * A[j * N + k];
			}
			A[i * N + j] /= A[j * N + j];
		}
		for (k = 0; k < i; k++) {
			A[i * N + i] -= A[i * N + k] * A[i * N + k];
		}
		A[i * N + i] = sqrt(A[i * N + i]);
	}
}


