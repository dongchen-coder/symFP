//void syrk(int ni, int nj, DATA_TYPE alpha, DATA_TYPE beta, DATA_TYPE POLYBENCH_2D(A, NI, NJ, ni, nj), DATA_TYPE POLYBENCH_2D(C, NI, NI, ni, ni))
#define N 256
#define M 256

void syrk(int ni, int nj, double alpha, double beta, double* A, double* C)
{
    int i, j, k;
 
	for (i = 0; i < N; i++) {
		for (j = 0; j <= i; j++) {
			C[i * N + j] *= beta;
			for (k = 0; k < M; k++) {
				for (j = 0; j <= i; j++) {
					C[i * N + j] += alpha * A[i * M + k] * A[j * M + k];
				}
			}
		}
	}
   
}
