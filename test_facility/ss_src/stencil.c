#ifndef DEBUG
#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET) && !defined(MEDIUM_DATASET)
    #define LARGE_DATASET
#endif
#ifdef MINI_DATASET
	#define NI 32
	#define NJ 32
#endif
#ifdef SMALL_DATASET
	#define NI 1024
	#define NJ 1024
#endif
#ifdef MEDIUM_DATASET
	#define NI 4096
	#define NJ 4096
#endif
#ifdef LARGE_DATASET
    #define NI 8192
    #define NJ 8192
#endif
#ifdef EXTRALARGE_DATASET
	#define NI 100000
	#define NJ 100000
#endif
#else
	#define NI 4
	#define NJ 4
#endif


int stencil(double* a, double* b) {

	for (int i = 0; i < NI; i++) {
		for (int j = 0; j < NJ; j++) {
			b[(NJ)*(i+1)+j+1] = a[(NJ)*(i+1)+j+1] 
						+ a[(NJ)*(i+1)+j+2] 
						+ a[(NJ)*(i+1)+j]
						+ a[(NJ)*i+j+1]
						+ a[(NJ)*(i+2)+j+1];
		}
	}

	return 0;
}


