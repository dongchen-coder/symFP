
#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif
# ifndef POLYBENCH_LOOP_ITERATION
#  define POLYBENCH_LOOP_ITERATION  1000
# endif
#ifdef MINI_DATASET
	#define M 16
	#define N 4
#endif
#ifdef SMALL_DATASET
	#define M 128
    #define N 128
#endif 
#ifdef MEDIUM_DATASET
	#define M 1024
	#define N 1024
#endif
#ifdef LARGE_DATASET
	#define M 8192
	#define N 8192
#endif
#ifdef EXTRALARGE_DATASET
	#define M 4096
	#define N 4096
#endif

#else
	#define M 4
	#define N 4
#endif

void trmm(double* A, double* B, double alpha) {

	int i, j, k;

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			for (k = i+1; k < M; k++) {
				B[i * N + j] += A[k * M + i] * B[k *N + j];
			}
			B[i * N + j] = alpha * B[i * N + j];
     	}
	}
}
