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

#ifndef BLOCKSIZE
	#define BLOCKSIZE NI / 4
#endif

void mm_tile(double * pA, double * pB, double * pC) {

	int i, j, k, ii, jj, kk;
	double sum;

	for (ii = 0; ii < NI; ii += BLOCKSIZE)
	{
		for( jj = 0; jj < NJ; jj += BLOCKSIZE )
		{
			for( kk = 0; kk < NK; kk += BLOCKSIZE )
			{
				for( i = ii; i < (ii + BLOCKSIZE); i += 1 )
				{
					for ( j = jj; j < (jj + BLOCKSIZE); j += 1)
					{
						sum = 0.0;
						for ( k = kk; k < (kk + BLOCKSIZE); k += 1)
						{
							sum += pA[NK*i+k] * pB[NJ*k+j];
						}
						pC[NJ*i+j] = sum;
					}
				}
			}
		}
	}
}