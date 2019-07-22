
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"

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
            rtTmpAccess(A_OFFSET + i * N + j, 0, 0);
            for (k = 0; k < j; k++) {
                w -= A[i * N + k] * A[k * N + j];
                rtTmpAccess(A_OFFSET + i * N + k, 1, 0);
                rtTmpAccess(A_OFFSET + k * N + j, 2, 0);
            }
            A[i * N + j] = w / A[j * N + j];
            rtTmpAccess(A_OFFSET + j * N + j, 3, 0);
            rtTmpAccess(A_OFFSET + i * N + j, 4, 0);
        }
        for (j = i; j < N; j++) {
            w = A[i * N + j];
            rtTmpAccess(A_OFFSET + i * N + j, 5, 0);
            for (k = 0; k < i; k++) {
                w -= A[i * N + k] * A[k * N + j];
                rtTmpAccess(A_OFFSET + i * N + k, 6, 0);
                rtTmpAccess(A_OFFSET + k * N + j, 7, 0);
            }
            A[i * N + j] = w;
            rtTmpAccess(A_OFFSET + i * N + j, 8, 0);
        }
    }
    
    for (i = 0; i < N; i++) {
        w = b[i];
        rtTmpAccess(B_OFFSET + i, 9, 1);
        for (j = 0; j < i; j++) {
            w -= A[i * N + j] * y[j];
            rtTmpAccess(A_OFFSET + i * N + j, 10, 0);
            rtTmpAccess(Y_OFFSET + j, 11, 2);
        }
        y[i] = w;
        rtTmpAccess(Y_OFFSET + i, 12, 2);
    }
    for (i = N-1; i >=0; i--) {
        w = y[i];
        rtTmpAccess(Y_OFFSET + i, 13, 2);
        for (j = i+1; j < N; j++) {
            w -= A[i * N + j] * x[j];
            rtTmpAccess(A_OFFSET + i * N + j, 14, 0);
            rtTmpAccess(X_OFFSET + j, 15, 3);
        }
        x[i] = w / A[i * N + i];
        rtTmpAccess(A_OFFSET + i * N + i, 16, 0);
        rtTmpAccess(X_OFFSET + i, 17, 3);
    }
}

int main() {

	double* A = (double *)malloc(N * N *sizeof(double));
	double* b = (double *)malloc(N * sizeof(double));
	double* y = (double *)malloc(N * sizeof(double));
	double* x = (double *)malloc(N * sizeof(double));

	ludcmp_trace(A, b, y, x);

	OSL_ref(0);

	return 0;
}
