
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"
#ifdef PAPI_TIMER
#   include "../utility/papi_timer.h"
#endif
#ifdef ORG
	#define M 1024
	#define N 1024
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define A_OFFSET 0
#define B_OFFSET M * M

void trmm_trace(double* A, double* B, double alpha) {
    
    int i, j, k;
    
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            for (k = i+1; k < M; k++) {
                B[i * N + j] += A[k * M + i] * B[k * N + j];
                rtTmpAccess(A_OFFSET + k * M + i, 0, 0);
                rtTmpAccess(B_OFFSET + k * N + j, 1, 1);
                rtTmpAccess(B_OFFSET + i * N + j, 2, 1);
                rtTmpAccess(B_OFFSET + i * N + j, 3, 1);
            }
            B[i * N + j] = alpha * B[i * N + j];
            rtTmpAccess(B_OFFSET + i * N + j, 4, 1);
            rtTmpAccess(B_OFFSET + i * N + j, 5, 1);
        }
    }
}

int main() {

	double* A = (double *)malloc(M * M * sizeof(double));
	double* B = (double *)malloc(M * N * sizeof(double));
	double alpha = 0.2;

#ifdef PAPI_TIMER
    PAPI_timer_init();
    PAPI_timer_start();
#endif
	trmm_trace(A, B, alpha);
    
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
