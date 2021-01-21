#ifndef DEBUG
#  if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#   define STANDARD_DATASET
#  endif
#  ifdef MINI_DATASET
#   define NQ 32
#   define NR 32
#   define NP 32
#  endif

#  ifdef SMALL_DATASET
#   define NQ 128
#   define NR 128
#   define NP 128
#  endif

#  ifdef STANDARD_DATASET /* Default if unspecified. */
#   define NQ 256
#   define NR 256
#   define NP 256
#  endif

#  ifdef LARGE_DATASET
#   define NQ 1024
#   define NR 1024
#   define NP 1024
#  endif

#  ifdef EXTRALARGE_DATASET
#   define NQ 4096
#   define NR 4096
#   define NP 4096
#  endif
#else
#	define NQ 4
#	define NR 4
#	define NP 4
#endif


void doitgen(double* sum, double* A, double* C4) {

	int r, q, p, s;

	for (r = 0; r < NR; r++) {
		for (q = 0; q < NQ; q++)  {
			for (p = 0; p < NP; p++)  {
				sum[r * NQ * NP + q * NP + p] = 0.0;
				for (s = 0; s < NP; s++) {
					sum[r * NQ * NP + q * NP + p] += A[r * NQ * NP + q * NP + s] * C4[s * NP + p];
				}
			}
			for (p = 0; p < NP; p++) {
				A[r * NQ * NP + q * NP + p] = sum[r * NQ * NP + q * NP + p];
			}
    	}
	}

}

