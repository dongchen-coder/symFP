

#define TSTEPS 10
#define N 1024

void jacobi_2d(double* A, double* B) {

	int t, i, j;

	for (t = 0; t < TSTEPS; t++) {
		for (i = 1; i < N - 1; i++)
			for (j = 1; j < N - 1; j++)
				B[i * N + j] = 0.2 * (A[i * N + j] + A[i * N + j-1] + A[i * N + 1+j] + A[(1+i) * N + j] + A[(i-1) * N + j]);
		for (i = 1; i < N - 1; i++)
			for (j = 1; j < N - 1; j++)
				A[i * N + j] = 0.2 * (B[i * N + j] + B[i * N + j-1] + B[i * N + 1+j] + B[(1+i) * N + j] + B[(i-1) * N + j]);
	}

}

