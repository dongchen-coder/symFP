#ifndef DEBUG
#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
    #define STANDARD_DATASET
#endif

/* Define the possible dataset sizes. */
#  ifdef MINI_DATASET
#   define N 32
#   define M 32
#  endif

#  ifdef SMALL_DATASET
#   define N 128
#   define M 128
#  endif

#  ifdef STANDARD_DATASET /* Default if unspecified. */
#   define N 1024
#   define M 1024
#  endif

#  ifdef LARGE_DATASET
#   define N 2048
#   define M 2048
#  endif

#  ifdef EXTRALARGE_DATASET
#   define N 4096
#   define M 4096
#  endif
#else
#	define N 4
# 	define M 4
#endif

void syr2k(double* A, double* B, double* C, double alpha, double beta) {

	int i, j, k;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			C[i* N + j] *= beta;
		}
    }
    for (i = 0; i < N; i++) {
    	for (j = 0; j < N; j++) {
    		for (k = 0; k < M; k++) {
    			C[i * N + j] += A[j * M + k] * alpha * B[i * M + k] + B[j * M + k] * alpha * A[i * N + k];
        	}
    	}
    }

}

