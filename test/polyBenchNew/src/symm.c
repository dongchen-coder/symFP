

#define M 1024
#define N 1024

void symm(double* A, double* B, double* C, double alpha, double beta) {
	
	int i, j, k;	
	double temp2;

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++ ) {
			temp2 = 0;
			for (k = 0; k < i; k++) {
				C[k * N + j] += alpha*B[i * N + j] * A[i * M + k];
				temp2 += B[k * N + j] * A[i * M + k];
			}
			C[i * N + j] = beta * C[i * N + j] + alpha*B[i * N + j] * A[i * M + i] + alpha * temp2;
		}
	}
	return;
}


