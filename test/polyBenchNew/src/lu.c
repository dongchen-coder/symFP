
#define N 1024

void lu(double* A) {

	int i, j, k;

	for (i = 0; i < N; i++) {
		for (j = 0; j <i; j++) {
			for (k = 0; k < j; k++) {
				A[i * N + j] -= A[i * N + k] * A[k * N + j];
			}
			A[i * N + j] /= A[j * N + j];
		}
		for (j = i; j < N; j++) {
			for (k = 0; k < i; k++) {
				A[i * N + j] -= A[i * N + k] * A[k * N + j];
			}
		}
	}
}

