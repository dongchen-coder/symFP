
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define N 1024
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define Y_OFFSET 0
#define R_OFFSET N
#define Z_OFFSET N + N

void durbin_trace(double *y, double* r, double* z) {
    
    int k, i;
    double alpha, beta, sum;
    
    /*
     y[0] = -r[0];
     beta = 1.0;
     alpha = -r[0];
     */
    for (k = 1; k < N; k++) {
        beta = (1-alpha*alpha)*beta;
        sum = 0.0;
        for (i=0; i<k; i++) {
            sum += r[k-i-1]*y[i];
            rtTmpAccess(R_OFFSET + k-i-1, 0, 0);
            rtTmpAccess(Y_OFFSET + i, 1, 1);
        }
        alpha = - (r[k] + sum)/beta;
        rtTmpAccess(R_OFFSET + k, 2, 0);
        
        for (i=0; i<k; i++) {
            z[i] = y[i] + alpha*y[k-i-1];
            rtTmpAccess(Y_OFFSET + i, 3, 1);
            rtTmpAccess(Y_OFFSET + k-i-1, 4, 1);
            rtTmpAccess(Z_OFFSET + i, 5, 2);
        }
        for (i=0; i<k; i++) {
            y[i] = z[i];
            rtTmpAccess(Z_OFFSET + i, 6, 2);
            rtTmpAccess(Y_OFFSET + i, 7, 1);
        }
        y[k] = alpha;
        rtTmpAccess(Y_OFFSET + k, 8, 1);
    }
    
}

int main() {

	double * y = (double *) malloc(N * sizeof(double));
	double * r = (double *) malloc(N * sizeof(double));
	double * z = (double *) malloc(N * sizeof(double));

	durbin_trace(y, r, z);

	OSL_ref(0);

	return 0;
}

