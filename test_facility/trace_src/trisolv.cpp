#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

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

#ifdef RD
    InitRD();
#endif
    
	trisolv_trace(x, b, L);

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


