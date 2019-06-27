#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

#ifdef ORG
	#define TSTEPS 10
	#define N 1024
#elif defined (TX)
#elif defined (FX)
#elif defined (EX)
#endif

#define A_OFFSET 0
#define B_OFFSET N * N

void jacobi_2d_trace(double* A, double* B) {

	int t, i, j;

	for (t = 0; t < TSTEPS; t++) {
		for (i = 1; i < N - 1; i++) {
			for (j = 1; j < N - 1; j++) {
				B[i * N + j] = 0.2 * (A[i * N + j] + A[i * N + j-1] + A[i * N + 1+j] + A[(1+i) * N + j] + A[(i-1) * N + j]);
				rtTmpAccess(A_OFFSET + i * N + j);
				rtTmpAccess(A_OFFSET + i * N + j-1);
				rtTmpAccess(A_OFFSET + i * N + 1+j);
				rtTmpAccess(A_OFFSET + (1+i) * N + j);
				rtTmpAccess(A_OFFSET + (i-1) * N + j);
				rtTmpAccess(B_OFFSET + i * N + j);
			}
		}
		for (i = 1; i < N - 1; i++) {
			for (j = 1; j < N - 1; j++) {
				A[i * N + j] = 0.2 * (B[i * N + j] + B[i * N + j-1] + B[i * N + 1+j] + B[(1+i) * N + j] + B[(i-1) * N + j]);
				rtTmpAccess(B_OFFSET + i * N + j);
				rtTmpAccess(B_OFFSET + i * N + j-1);
				rtTmpAccess(B_OFFSET + i * N + 1+j);
				rtTmpAccess(B_OFFSET + (1+i) * N + j);
				rtTmpAccess(B_OFFSET + (i-1) * N + j);
				rtTmpAccess(A_OFFSET + i * N + j);
			}
		}
	}
}

int main() {

	double* A = (double *)malloc(N * N * sizeof(double));
	double* B = (double *)malloc(N * N * sizeof(double));

#ifdef RD
    InitRD();
#endif
    
	jacobi_2d_trace(A, B);

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
