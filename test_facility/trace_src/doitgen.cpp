#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

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

#define SUM_OFFSET 0
#define A_OFFSET NP
#define C4_OFFSET NP + NR * NQ * NP

void doitgen_trace(double* sum, double* A, double* C4) {

	int r, q, p, s;

	for (r = 0; r < NR; r++) {
		for (q = 0; q < NQ; q++)  {
			for (p = 0; p < NP; p++)  {
				sum[p] = 0.0;
				rtTmpAccess(SUM_OFFSET + p);
				for (s = 0; s < NP; s++) {
					sum[p] += A[r * NQ * NP + q * NP + s] * C4[s * NP + p];
					rtTmpAccess(A_OFFSET + r * NQ * NP + q * NP + s);
					rtTmpAccess(C4_OFFSET + s * NP + p);
					rtTmpAccess(SUM_OFFSET + p);
					rtTmpAccess(SUM_OFFSET + p);
				}
			}
			for (p = 0; p < NP; p++) {
				A[r * NQ * NP + q * NP + p] = sum[p];
				rtTmpAccess(SUM_OFFSET + p);
				rtTmpAccess(A_OFFSET + r * NQ * NP + q * NP + p);
			}
    	}
	}

}

int main() {
	
	double * sum = (double *) malloc(NP * sizeof(double));
	double * A = (double *) malloc(NR * NQ * NP * sizeof(double));
	double * C4 = (double *) malloc(NP * NP);
	
#ifdef RD
    InitRD();
#endif
    
	doitgen_trace(sum, A, C4);

#ifdef PROFILE_RT
    dumpRtTmp();
    RTtoMR_AET();
    dumpMR();
#endif
    
#ifdef RD
    FiniRD();
#endif
    free(sum);
    free(A);
    free(C4);

	return 0;
}



