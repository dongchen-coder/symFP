
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define N 1024
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define A_OFFSET 0

void lu_trace(double* A) {
    
    int i, j, k;
    
    for (i = 0; i < N; i++) {
        for (j = 0; j <i; j++) {
            for (k = 0; k < j; k++) {
                A[i * N + j] -= A[i * N + k] * A[k * N + j];
                rtTmpAccess(A_OFFSET + i * N + k, 0, 0);
                rtTmpAccess(A_OFFSET + k * N + j, 1, 0);
                rtTmpAccess(A_OFFSET + i * N + j, 2, 0);
                rtTmpAccess(A_OFFSET + i * N + j, 3, 0);
            }
            A[i * N + j] /= A[j * N + j];
            rtTmpAccess(A_OFFSET + j * N + j, 4, 0);
            rtTmpAccess(A_OFFSET + i * N + j, 5, 0);
            rtTmpAccess(A_OFFSET + i * N + j, 6, 0);
        }
        for (j = i; j < N; j++) {
            for (k = 0; k < i; k++) {
                A[i * N + j] -= A[i * N + k] * A[k * N + j];
                rtTmpAccess(A_OFFSET + i * N + k, 7, 0);
                rtTmpAccess(A_OFFSET + k * N + j, 8, 0);
                rtTmpAccess(A_OFFSET + i * N + j, 9, 0);
                rtTmpAccess(A_OFFSET + i * N + j, 10, 0);
            }
        }
    }
}

int main() {

	double* A = (double *)malloc(N * N * sizeof(double));

	lu_trace(A);

	OSL_ref(0);

	return 0;
}
