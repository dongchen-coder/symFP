
#define TSTEPS 10
#define N 1024

void jacobi_1d(double* A, double* B) {

	int t, i;

	for (i = 1; i < N - 1; i++)
		B[i] = 0.33333 * (A[i-1] + A[i] + A[i + 1]);
	for (i = 1; i < N - 1; i++)
		A[i] = 0.33333 * (B[i-1] + B[i] + B[i + 1]);

}


