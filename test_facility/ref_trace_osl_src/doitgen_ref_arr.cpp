
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define NR 64
	#define NQ 64
	#define NP 64
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define SUM_OFFSET 0
#define A_OFFSET NP
#define C4_OFFSET NP + NR * NQ * NP

void doitgen_trace(double* sum, double* A, double* C4) {
    
    int r, q, p, s;
    
    for (r = 0; r < NR; r++) {
        for (q = 0; q < NQ; q++)  {
            for (p = 0; p < NP; p++)  {
                sum[p] = 0.0;
                rtTmpAccess(SUM_OFFSET + p, 0, 0);
                for (s = 0; s < NP; s++) {
                    sum[p] += A[r * NQ * NP + q * NP + s] * C4[s * NP + p];
                    rtTmpAccess(A_OFFSET + r * NQ * NP + q * NP + s, 1, 1);
                    rtTmpAccess(C4_OFFSET + s * NP + p, 2, 2);
                    rtTmpAccess(SUM_OFFSET + p, 3, 0);
                    rtTmpAccess(SUM_OFFSET + p, 4, 0);
                }
            }
            for (p = 0; p < NP; p++) {
                A[r * NQ * NP + q * NP + p] = sum[p];
                rtTmpAccess(SUM_OFFSET + p, 5, 0);
                rtTmpAccess(A_OFFSET + r * NQ * NP + q * NP + p, 6, 1);
            }
        }
    }
    
}

int main() {
	
	double * sum = (double *) malloc(NP * sizeof(double));
	double * A = (double *) malloc(NR * NQ * NP * sizeof(double));
	double * C4 = (double *) malloc(NP * NP);
	
	doitgen_trace(sum, A, C4);

	OSL_ref(0);

	return 0;
}


