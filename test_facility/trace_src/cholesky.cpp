#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif


#include <math.h>

#ifdef ORG
	#define N 1024
#elif defined (TX)
#elif defined (FX)
#elif defined (EX)
#endif

#define A_OFFSET 0

void cholesky_trace(double* A) {

	int i, j, k;

	for (i = 0; i < N; i++) {
		for (j = 0; j < i; j++) {
			for (k = 0; k < j; k++) {
				A[i * N + j] -= A[i * N + k] * A[j * N + k];
				rtTmpAccess(A_OFFSET + i * N + k);
				rtTmpAccess(A_OFFSET + j * N + k);
				rtTmpAccess(A_OFFSET + i * N + j);
				rtTmpAccess(A_OFFSET + i * N + j);
			}
			A[i * N + j] /= A[j * N + j];
			rtTmpAccess(A_OFFSET + j * N + j);
			rtTmpAccess(A_OFFSET + i * N + j);
			rtTmpAccess(A_OFFSET + i * N + j);
		}
		for (k = 0; k < i; k++) {
			A[i * N + i] -= A[i * N + k] * A[i * N + k];
			rtTmpAccess(A_OFFSET + i * N + k);
			rtTmpAccess(A_OFFSET + i * N + k);
			rtTmpAccess(A_OFFSET + i * N + i);
			rtTmpAccess(A_OFFSET + i * N + i);
		}
		A[i * N + i] = sqrt(A[i * N + i]);
		rtTmpAccess(A_OFFSET + i * N + i);
		rtTmpAccess(A_OFFSET + i * N + i);
	}
}

int main() {

	double * A = (double *) malloc(N * N * sizeof(double));
	
#ifdef RD
    InitRD();
#endif
    
	cholesky_trace(A);
	
#ifdef PROFILE_RT
    dumpRtTmp();
    RTtoMR_AET();
    dumpMR();
#endif
    
#ifdef RD
    FiniRD();
#endif

	return 0;
}


