#include "../utility/rt.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define N 1024
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define X_OFFSET 0
#define B_OFFSET N
#define L_OFFSET N + N

void trisolv_trace(double* x, double* b, double* L) {

	int i, j;

	for (i = 0; i < N; i++) {
		x[i] = b[i];
		rtTmpAccess(X_OFFSET + i);
		rtTmpAccess(B_OFFSET + i);
		for (j = 0; j <i; j++) {
			x[i] -= L[i * N + j] * x[j];
			rtTmpAccess(L_OFFSET + i * N + j);
			rtTmpAccess(X_OFFSET + j);
			rtTmpAccess(X_OFFSET + i);
			rtTmpAccess(X_OFFSET + i);
		}
		x[i] = x[i] / L[i * N + i];
		rtTmpAccess(X_OFFSET + i);
		rtTmpAccess(L_OFFSET + i * N + i);
		rtTmpAccess(X_OFFSET + i);
		rtTmpAccess(X_OFFSET + i);
	}

}

int main() {

	double* x = (double *)malloc(N * sizeof(double));
	double* b = (double *)malloc(N * sizeof(double));
	double* L = (double *)malloc(N * N * sizeof(double));

	trisolv_trace(x, b, L);

	dumpRtTmp();
    RTtoMR_AET();
    dumpMR();

	return 0;
}


