#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif
# ifndef POLYBENCH_LOOP_ITERATION
#  define POLYBENCH_LOOP_ITERATION  1000
# endif
#ifdef MINI_DATASET
	#define NI 32
	#define NJ 32
	#define NK 32
	#define NL 32
#endif
#ifdef SMALL_DATASET
	#define NI 128
    #define NJ 128
    #define NL 128
    #define NK 128
#endif 
#ifdef MEDIUM_DATASET
	#define NI 1024
	#define NJ 1024
	#define NK 1024
	#define NL 1024
#endif
#ifdef LARGE_DATASET
	#define NI 8192
	#define NJ 8192
	#define NL 8192
	#define NK 8192
#endif
#ifdef EXTRALARGE_DATASET
	#define NI 4096
	#define NJ 4096
	#define NL 4096
	#define NK 4096
#endif

#else
	#define NI 4
	#define NJ 4
	#define NK 4
	#define NL 4
#endif

void mm2(double* tmp, double* A, double* B, double* C, double* D, double alpha, double beta) {

	int i, j, k;

	for (i = 0; i < NI; i++) {
		for (j = 0; j < NJ; j++) {
			tmp[i * NJ + j] = 0.0;
			for (k = 0; k < NK; ++k) {
				tmp[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
      		}
		}
	}
	for (i = 0; i < NI; i++) {
		for (j = 0; j < NL; j++) {
			D[i * NL + j] *= beta;
			for (k = 0; k < NJ; ++k) {
				D[i * NL + j] += tmp[i * NJ + k] * C[k * NL + j];
      		}
		}
	}
}
