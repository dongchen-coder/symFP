
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"
#ifdef PAPI_TIMER
#   include "../utility/papi_timer.h"
#endif
#include <math.h>

#ifdef ORG
	#define N 1024
#elif defined (TX)
#elif defined (FX)
#elif defined (EX)
#endif

#define A_OFFSET 0

void cholesky_trace(double* A) {
    
    int i, j, k;
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < i; j++) {
            for (k = 0; k < j; k++) {
                A[i * N + j] -= A[i * N + k] * A[j * N + k];
                rtTmpAccess(A_OFFSET + i * N + k, 0, 0);
                rtTmpAccess(A_OFFSET + j * N + k, 1, 0);
                rtTmpAccess(A_OFFSET + i * N + j, 2, 0);
                rtTmpAccess(A_OFFSET + i * N + j, 3, 0);
            }
            A[i * N + j] /= A[j * N + j];
            rtTmpAccess(A_OFFSET + j * N + j, 4, 0);
            rtTmpAccess(A_OFFSET + i * N + j, 5, 0);
            rtTmpAccess(A_OFFSET + i * N + j, 6, 0);
        }
        for (k = 0; k < i; k++) {
            A[i * N + i] -= A[i * N + k] * A[i * N + k];
            rtTmpAccess(A_OFFSET + i * N + k, 7, 0);
            rtTmpAccess(A_OFFSET + i * N + k, 8, 0);
            rtTmpAccess(A_OFFSET + i * N + i, 9, 0);
            rtTmpAccess(A_OFFSET + i * N + i, 10, 0);
        }
        A[i * N + i] = sqrt(A[i * N + i]);
        rtTmpAccess(A_OFFSET + i * N + i, 11, 0);
        rtTmpAccess(A_OFFSET + i * N + i, 12, 0);
    }
}

int main() {

	double * A = (double *) malloc(N * N * sizeof(double));
	
#ifdef PAPI_TIMER
    PAPI_timer_init();
    PAPI_timer_start();
#endif
	cholesky_trace(A);
#ifdef PAPI_TIMER
    PAPI_timer_end();
    PAPI_timer_print();
#endif

    printf("CARL Lease Assignment:\n");
#ifdef PAPI_TIMER
    PAPI_timer_start();
#endif
    OSL_ref(0);
#ifdef PAPI_TIMER
    PAPI_timer_end();
    PAPI_timer_print();
#endif
	return 0;
}

