
#define TSTEPS 10
#define N 1024

void seidel_2d(double* A) {

	int t, i, j;

	for (t = 0; t <= TSTEPS - 1; t++)
		for (i = 1; i<= N - 2; i++)
			for (j = 1; j <= N - 2; j++)
				A[i * N + j] = (A[(i-1) * N + j-1] + A[(i-1) * N + j] + A[(i-1) * N + j+1] + A[i * N + j-1] + A[i * N + j] + A[i * N + j+1] + A[(i+1) * N + j-1] + A[(i+1) * N + j] + A[(i+1) * N + j+1]) / 9.0;

}
