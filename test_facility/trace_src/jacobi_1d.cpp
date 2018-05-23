#include "../utility/rt.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define TSTEPS 10
	#define N 1024
#elif defined (TX)
#elif defined (FX)
#elif defined (EX)
#endif

#define A_OFFSET 0
#define B_OFFSET N

void jacobi_1d_trace(double* A, double* B) {

	int t, i;

	for (t = 0; t < TSTEPS; t++) {
		for (i = 1; i < N - 1; i++) {
			B[i] = 0.33333 * (A[i-1] + A[i] + A[i + 1]);
			rtTmpAccess(A_OFFSET + i - 1);
			rtTmpAccess(A_OFFSET + i);
			rtTmpAccess(A_OFFSET + i + 1);
			rtTmpAccess(B_OFFSET + i);
		}
		for (i = 1; i < N - 1; i++) {
			A[i] = 0.33333 * (B[i-1] + B[i] + B[i + 1]);
			rtTmpAccess(B_OFFSET + i - 1);
			rtTmpAccess(B_OFFSET + i);
			rtTmpAccess(B_OFFSET + i + 1);
			rtTmpAccess(A_OFFSET + i);
		}
	}
}

int main() {

	double* A = (double *)malloc(N * sizeof(double));
	double* B = (double *)malloc(N * sizeof(double));

	jacobi_1d_trace(A, B);

	dumpRtTmp();
    RTtoMR_AET();
    dumpMR();
	
	return 0;
}

