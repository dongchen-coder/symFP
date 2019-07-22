
#include "../utility/reference_lease.h"
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
        rtTmpAccess(B_OFFSET + i, 0, 1);
        rtTmpAccess(X_OFFSET + i, 1, 0);
        
        for (j = 0; j <i; j++) {
            x[i] -= L[i * N + j] * x[j];
            rtTmpAccess(L_OFFSET + i * N + j, 2, 2);
            rtTmpAccess(X_OFFSET + j, 3, 0);
            rtTmpAccess(X_OFFSET + i, 4, 0);
            rtTmpAccess(X_OFFSET + i, 5, 0);
        }
        x[i] = x[i] / L[i * N + i];
        rtTmpAccess(X_OFFSET + i, 6, 0);
        rtTmpAccess(L_OFFSET + i * N + i, 7, 2);
        rtTmpAccess(X_OFFSET + i, 8, 0);
    }
    
}

int main() {

	double* x = (double *)malloc(N * sizeof(double));
	double* b = (double *)malloc(N * sizeof(double));
	double* L = (double *)malloc(N * N * sizeof(double));

	trisolv_trace(x, b, L);

	OSL_ref(0);

	return 0;
}

