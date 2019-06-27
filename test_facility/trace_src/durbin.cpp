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
			rtTmpAccess(R_OFFSET + k-i-1);
			rtTmpAccess(Y_OFFSET + i);
		}
		alpha = - (r[k] + sum)/beta;
		rtTmpAccess(R_OFFSET + k);

		for (i=0; i<k; i++) {
			z[i] = y[i] + alpha*y[k-i-1];
			rtTmpAccess(Y_OFFSET + i);
			rtTmpAccess(Y_OFFSET + k-i-1);
			rtTmpAccess(Z_OFFSET + i);
		}
		for (i=0; i<k; i++) {
			y[i] = z[i];
			rtTmpAccess(Z_OFFSET + i);
			rtTmpAccess(Y_OFFSET + i);
		}
		y[k] = alpha;
		rtTmpAccess(Y_OFFSET + k);
	}	

}

int main() {

	double * y = (double *) malloc(N * sizeof(double));
	double * r = (double *) malloc(N * sizeof(double));
	double * z = (double *) malloc(N * sizeof(double));

#ifdef RD
    InitRD();
#endif
    
	durbin_trace(y, r, z);

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


