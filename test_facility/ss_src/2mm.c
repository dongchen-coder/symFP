#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define STANDARD_DATASET
# endif
#ifdef MINI_DATASET
	#define NI 32
	#define NJ 32
	#define NK 32
	#define NL 32
#endif
#ifdef SMALL_DATASET
	#define NI 256
    #define NJ 256
    #define NL 256
    #define NK 256
#endif 
#ifdef STANDARD_DATASET
	#define NI 1024
	#define NJ 1024
	#define NK 1024
	#define NL 1024
#endif
#ifdef LARGE_DATASET
	#define NI 2048
	#define NJ 2048
	#define NL 2048
	#define NK 2048
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
