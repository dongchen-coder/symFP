
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"


#ifdef ORG
	#define TSTEPS 10
	#define N 256
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define A_OFFSET 0
#define B_OFFSET N * N * N

void heat_3d_trace(double* B, double* A) {
    
    int t, i, j, k;
    
    for (t = 1; t <= TSTEPS; t++) {
        for (i = 1; i < N-1; i++) {
            for (j = 1; j < N-1; j++) {
                for (k = 1; k < N-1; k++) {
                    B[i * N * N + j * N + k] =   0.125 * (A[(i+1) * N * N + j * N + k] - 2.0 * A[i * N * N + j * N + k] + A[(i-1) * N * N + j * N + k])
                    + 0.125 * (A[i * N * N + (j+1) * N + k] - 2.0 * A[i * N * N + j * N + k] + A[i * N * N + (j-1) * N + k])
                    + 0.125 * (A[i * N * N + j * N + k+1] - 2.0 * A[i * N * N + j * N + k] + A[i * N * N + j * N + k-1])
                    + A[i * N * N + j * N + k];
                    rtTmpAccess(A_OFFSET + (i+1) * N * N + j * N + k, 0, 0);
                    rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 1, 0);
                    rtTmpAccess(A_OFFSET + (i-1) * N * N + j * N + k, 2, 0);
                    rtTmpAccess(A_OFFSET + i * N * N + (j+1) * N + k, 3, 0);
                    rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 4, 0);
                    rtTmpAccess(A_OFFSET + i * N * N + (j-1) * N + k, 5, 0);
                    rtTmpAccess(A_OFFSET + i * N * N + j * N + k+1, 6, 0);
                    rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 7, 0);
                    rtTmpAccess(A_OFFSET + i * N * N + j * N + k-1, 8, 0);
                    rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 9, 0);
                    rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 10, 1);
                }
            }
        }
        for (i = 1; i < N-1; i++) {
            for (j = 1; j < N-1; j++) {
                for (k = 1; k < N-1; k++) {
                    A[i * N * N + j * N + k] =   0.125 * (B[(i+1) * N * N + j * N + k] - 2.0 * B[i * N * N + j * N + k] + B[(i-1) * N * N + j * N + k])
                    + 0.125 * (B[i * N * N + (j+1) * N + k] - 2.0 * B[i * N * N + j * N + k] + B[i * N * N + (j-1) * N + k])
                    + 0.125 * (B[i * N * N + j * N + k+1] - 2.0 * B[i * N * N + j * N + k] + B[i * N * N + j * N + k-1])
                    + B[i * N * N + j * N + k];
                    rtTmpAccess(B_OFFSET + (i+1) * N * N + j * N + k, 11, 1);
                    rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 12, 1);
                    rtTmpAccess(B_OFFSET + (i-1) * N * N + j * N + k, 13, 1);
                    rtTmpAccess(B_OFFSET + i * N * N + (j+1) * N + k, 14, 1);
                    rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 15, 1);
                    rtTmpAccess(B_OFFSET + i * N * N + (j-1) * N + k, 16, 1);
                    rtTmpAccess(B_OFFSET + i * N * N + j * N + k+1, 17, 1);
                    rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 18, 1);
                    rtTmpAccess(B_OFFSET + i * N * N + j * N + k-1, 19, 1);
                    rtTmpAccess(B_OFFSET + i * N * N + j * N + k, 20, 1);
                    rtTmpAccess(A_OFFSET + i * N * N + j * N + k, 21, 0);
                }
            }
        }
    }
}

int main() {

	double* A = (double *)malloc(N * N * N * sizeof(double));
	double* B = (double *)malloc(N * N * N * sizeof(double));

	heat_3d_trace(B, A);

	OSL_ref(0);

	return 0;
}
