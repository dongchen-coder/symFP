

#define TSTEPS 10
#define N 256

void heat_3d(double* B, double* A) {

	int t, i, j, k;

	for (t = 1; t <= TSTEPS; t++) {
		for (i = 1; i < N-1; i++) {
			for (j = 1; j < N-1; j++) {
				for (k = 1; k < N-1; k++) {
					B[i * N * N + j * N + k] =   0.125 * (A[(i+1) * N * N + j * N + k] - 2.0 * A[i * N * N + j * N + k] + A[(i-1) * N * N + j * N + k])
                                 + 0.125 * (A[i * N * N + (j+1) * N + k] - 2.0 * A[i * N * N + j * N + k] + A[i * N * N + (j-1) * N + k])
                                 + 0.125 * (A[i * N * N + j * N + k+1] - 2.0 * A[i * N * N + j * N + k] + A[i * N * N + j * N + k-1])
                                 + A[i * N * N + j * N + k];
                }
            }
        }
        for (i = 1; i < N-1; i++) {
           for (j = 1; j < N-1; j++) {
               for (k = 1; k < N-1; k++) {
                   A[i * N * N + j * N + k] =   0.125 * (B[(i+1) * N * N + j * N + k] - 2.0 * B[i * N * N + j * N + k] + B[(i-1) * N * N + j * N + k])
                                + 0.125 * (B[i * N * N + (j+1) * N + k] - 2.0 * B[i * N * N + j * N + k] + B[i * N * N + (j-1) * N + k])
                                + 0.125 * (B[i * N * N + j * N + k+1] - 2.0 * B[i * N * N + j * N + k] + B[i * N * N + j * N + k-1])
                                + B[i * N * N + j * N + k];
               }
           }
       }
    }

}
