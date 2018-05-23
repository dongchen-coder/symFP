#include "../utility/rt.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define TSTEPS 10
	#define N 1024
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define A_OFFSET 0

void seidel_2d_trace(double* A) {

	int t, i, j;

	for (t = 0; t <= TSTEPS - 1; t++)
		for (i = 1; i<= N - 2; i++)
			for (j = 1; j <= N - 2; j++) {
				A[i * N + j] = (A[(i-1) * N + j-1] + A[(i-1) * N + j] + A[(i-1) * N + j+1] + A[i * N + j-1] + A[i * N + j] + A[i * N + j+1] + A[(i+1) * N + j-1] + A[(i+1) * N + j] + A[(i+1) * N + j+1]) / 9.0;
				rtTmpAccess(A_OFFSET + (i-1) * N + j-1);
				rtTmpAccess(A_OFFSET + (i-1) * N + j);
				rtTmpAccess(A_OFFSET + (i-1) * N + j+1);
				rtTmpAccess(A_OFFSET + i * N + j-1);
				rtTmpAccess(A_OFFSET + i * N + j);
				rtTmpAccess(A_OFFSET + i * N + j+1);
				rtTmpAccess(A_OFFSET + (i+1) * N + j-1);
				rtTmpAccess(A_OFFSET + (i+1) * N + j);
				rtTmpAccess(A_OFFSET + (i+1) * N + j+1);
				rtTmpAccess(A_OFFSET + i * N + j);
			}
}

int main() {

	double* A = (double *)malloc(N * N * sizeof(double));

	seidel_2d_trace(A);

	dumpRtTmp();
    RTtoMR_AET();
    dumpMR();

	return 0;
}

