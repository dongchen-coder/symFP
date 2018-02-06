
#define N 1024

void ludcmp(double* A, double* b, double* y, double* x) {

	int i, j, k;

	double w;

	for (i = 0; i < N; i++) {
		for (j = 0; j <i; j++) {
			w = A[i * N + j];
			for (k = 0; k < j; k++) {
				w -= A[i * N + k] * A[k * N + j];
			}
			A[i * N + j] = w / A[j * N + j];
		}
		for (j = i; j < N; j++) {
			w = A[i * N + j];
			for (k = 0; k < i; k++) {
				w -= A[i * N + k] * A[k * N + j];
			}
			A[i * N + j] = w;
		}
	}

	for (i = 0; i < N; i++) {
		w = b[i];
		for (j = 0; j < i; j++)
			w -= A[i * N + j] * y[j];
		y[i] = w;
	}
	for (i = N-1; i >=0; i--) {
		w = y[i];
		for (j = i+1; j < N; j++)
			w -= A[i * N + j] * x[j];
		x[i] = w / A[i * N + i];
	}

}


