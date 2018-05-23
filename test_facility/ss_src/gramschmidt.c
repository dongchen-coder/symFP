
#define N 1024
#define M 1024

void gramschmidt(double* A, double* R, double* Q) {

	int k, i, j;
	double nrm;

	for (k = 0; k < N; k++) {
		nrm = 0.0;
		for (i = 0; i < M; i++) {
			nrm += A[i * N + k] * A[i * N + k];
		}
		R[k * N + k] = sqrt(nrm);
		for (i = 0; i < M; i++) {
			Q[i * N + k] = A[i * N + k] / R[k * N + k];
		}
		for (j = k + 1; j < N; j++) {
			R[k * N + j] = 0.0;
			for (i = 0; i < M; i++) {
				R[k * N + j] += Q[i * N + k] * A[i * N + j];
			}
			for (i = 0; i < M; i++) {
        		A[i * N + j] = A[i * N + j] - Q[i * N + k] * R[k * N + j];
			}
		}
    }
}

