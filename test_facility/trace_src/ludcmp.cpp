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

#define A_OFFSET 0
#define B_OFFSET N * N
#define Y_OFFSET N * N + N
#define X_OFFSET N * N + N + N

void ludcmp_trace(double* A, double* b, double* y, double* x) {

	int i, j, k;

	double w;

	for (i = 0; i < N; i++) {
		for (j = 0; j <i; j++) {
			w = A[i * N + j];
			rtTmpAccess(A_OFFSET + i * N + j);
			for (k = 0; k < j; k++) {
				w -= A[i * N + k] * A[k * N + j];
				rtTmpAccess(A_OFFSET + i * N + k);
				rtTmpAccess(A_OFFSET + k * N + j);
			}
			A[i * N + j] = w / A[j * N + j];
			rtTmpAccess(A_OFFSET + j * N + j);
			rtTmpAccess(A_OFFSET + i * N + j);
		}
		for (j = i; j < N; j++) {
			w = A[i * N + j];
			rtTmpAccess(A_OFFSET + i * N + j);
			for (k = 0; k < i; k++) {
				w -= A[i * N + k] * A[k * N + j];
				rtTmpAccess(A_OFFSET + i * N + k);
				rtTmpAccess(A_OFFSET + k * N + j);
			}
			A[i * N + j] = w;
			rtTmpAccess(A_OFFSET + i * N + j);
		}
	}

	for (i = 0; i < N; i++) {
		w = b[i];
		rtTmpAccess(B_OFFSET + i);
		for (j = 0; j < i; j++) {
			w -= A[i * N + j] * y[j];
			rtTmpAccess(A_OFFSET + i * N + j);
			rtTmpAccess(Y_OFFSET + j);
		}
		y[i] = w;
		rtTmpAccess(Y_OFFSET + i);
	}
	for (i = N-1; i >=0; i--) {
		w = y[i];
		rtTmpAccess(Y_OFFSET + i);
		for (j = i+1; j < N; j++) {
			w -= A[i * N + j] * x[j];
			rtTmpAccess(A_OFFSET + i * N + j);
			rtTmpAccess(X_OFFSET + j);
		}
		x[i] = w / A[i * N + i];
		rtTmpAccess(A_OFFSET + i * N + i);
		rtTmpAccess(X_OFFSET + i);
	}
}

int main() {

	double* A = (double *)malloc(N * N *sizeof(double));
	double* b = (double *)malloc(N * sizeof(double));
	double* y = (double *)malloc(N * sizeof(double));
	double* x = (double *)malloc(N * sizeof(double));

#ifdef RD
    InitRD();
#endif
    
	ludcmp_trace(A, b, y, x);

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
