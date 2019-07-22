
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"

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
                rtTmpAccess(A_OFFSET + i * N + j, 0, 0);
                rtTmpAccess(A_OFFSET + i * N + j-1, 1, 0);
                rtTmpAccess(A_OFFSET + i * N + 1+j, 2, 0);
                rtTmpAccess(A_OFFSET + (1+i) * N + j, 3, 0);
                rtTmpAccess(A_OFFSET + (i-1) * N + j, 4, 0);
                rtTmpAccess(B_OFFSET + i * N + j, 5, 1);
            }
        }
        for (i = 1; i < N - 1; i++) {
            for (j = 1; j < N - 1; j++) {
                A[i * N + j] = 0.2 * (B[i * N + j] + B[i * N + j-1] + B[i * N + 1+j] + B[(1+i) * N + j] + B[(i-1) * N + j]);
                rtTmpAccess(B_OFFSET + i * N + j, 6, 1);
                rtTmpAccess(B_OFFSET + i * N + j-1, 7, 1);
                rtTmpAccess(B_OFFSET + i * N + 1+j, 8, 1);
                rtTmpAccess(B_OFFSET + (1+i) * N + j, 9, 1);
                rtTmpAccess(B_OFFSET + (i-1) * N + j, 10, 1);
                rtTmpAccess(A_OFFSET + i * N + j, 11, 0);
            }
        }
    }
}

int main() {

	double* A = (double *)malloc(N * N * sizeof(double));
	double* B = (double *)malloc(N * N * sizeof(double));

	jacobi_2d_trace(A, B);

	OSL_ref(0);

	return 0;
}
