//void syrk(int ni, int nj, DATA_TYPE alpha, DATA_TYPE beta, DATA_TYPE POLYBENCH_2D(A, NI, NJ, ni, nj), DATA_TYPE POLYBENCH_2D(C, NI, NI, ni, ni))
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


void syrk(double alpha, double beta, double* A, double* C)
{
    int i, j, k;

    for (i = 0; i < N; i++) {
    	for (j = 0; j < N; j++) {
    		C[i * N + j] *= beta;
    	}
	}

    for (i = 0; i < N; i++) {
    	for (j = 0; j < N; j++) {
    		for (k = 0; k < M; k++) {
    			C[i * N + j] += alpha * A[i * M + k] * A[j * M + k];
    		}
    	}
    }

   
}
