#ifndef DEBUG
#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
    #define STANDARD_DATASET
#endif
#ifdef MINI_DATASET
	#define NI 64
	#define NJ 64
#endif
#ifdef SMALL_DATASET
	#define NI 1024
	#define NJ 1024
#endif
#ifdef STANDARD_DATASET
	#define NI 2048
	#define NJ 2048
#endif
#ifdef EXTRALARGE_DATASET
    #define NI 4096
    #define NJ 4096
#endif
#ifdef EXTRALARGE_DATASET
	#define NI 8192
	#define NJ 8192
#endif
#else
	#define NI 4
	#define NJ 4
#endif


int stencil(double* a, double* b) {

	for (int i = 1; i < NI + 1; i++) {
		for (int j = 1; j < NJ + 1; j++) {
			b[(NJ+2)*i+j] = a[(NJ+2)*i+j] 
						+ a[(NJ+2)*i+j+1] 
						+ a[(NJ+2)*i+j-1]
						+ a[(NJ+2)*(i-1)+j]
						+ a[(NJ+2)*(i+1)+j];
		}
	}

	return 0;
}


