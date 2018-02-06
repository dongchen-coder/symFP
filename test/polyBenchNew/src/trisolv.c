

#define N 1024

void trisolv(double* x, double* b, double* L) {

	int i, j;

	for (i = 0; i < N; i++) {
		x[i] = b[i];
		for (j = 0; j <i; j++)
			x[i] -= L[i * N + j] * x[j];
		x[i] = x[i] / L[i * N + i];
	}

}
