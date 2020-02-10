#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define STANDARD_DATASET
# endif
#ifdef MINI_DATASET
	#define NI 32
	#define NJ 32
	#define NK 32
#endif
#ifdef SMALL_DATASET
	#define NI 128
    #define NJ 128
    #define NK 128
#endif 
#ifdef STANDARD_DATASET
	#define NI 1024
	#define NJ 1024
	#define NK 1024
#endif
#ifdef LARGE_DATASET
	#define NI 2048
	#define NJ 2048
	#define NK 2048
#endif
#ifdef EXTRALARGE_DATASET
	#define NI 4096
	#define NJ 4096
	#define NK 4096
#endif

#else
	#define NI 4
	#define NJ 4
	#define NK 4
#endif

void mm(double* pA, double* pB, double* pC) {

	int i, j, k;

	for( i = 0; i < NI; i += 1 ) {
		for( j = 0; j < NJ; j += 1 ) {
            pC[NJ*i+j] = 0.0;
			for( k = 0; k < NK; k += 1 ) {
				pC[NJ*i+j] += pA[NK*i+k] * pB[NJ*k+j];
			}
		}
	}
}
