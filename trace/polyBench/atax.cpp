#include "../utility/rt.h"

#define NX 1024
#define NY 1024

#define A_OFFSET 0
#define X_OFFSET NX * NY
#define Y_OFFSET NX * NY + NY
#define TMP_OFFSET NX * NY + NY + NY


void atax_cpu_trace(double* A, double* x, double* y, double* tmp, unsigned int dim_size_x, unsigned int dim_size_y) {
    
	int i,j;

    for (i= 0; i < NY; i++)
    {
        y[i] = 0;
    	rtTmpAccess(Y_OFFSET + i);
	}

    for (i = 0; i < NX; i++)
    {
            tmp[i] = 0;
			rtTmpAccess(TMP_OFFSET + i);

            for (j = 0; j < NY; j++)
            {
                tmp[i] = tmp[i] + A[i * NY + j] * x[j];
            	rtTmpAccess(TMP_OFFSET + i);
				rtTmpAccess(A_OFFSET + i * NY + j);
				rtTmpAccess(X_OFFSET + j);
				rtTmpAccess(TMP_OFFSET + i);
			}

            for (j = 0; j < NY; j++)
            {
                y[j] = y[j] + A[i * NY + j] * tmp[i];
            	rtTmpAccess(Y_OFFSET + j);
				rtTmpAccess(A_OFFSET + i * NY + j);
				rtTmpAccess(TMP_OFFSET + i);
				rtTmpAccess(Y_OFFSET + j);
			}
    }


	return;
}



int main() {

    double* A = (double*)malloc( (NX * NY) * sizeof(double) );
    double* x = (double*)malloc( NX * sizeof(double) );
    double* y = (double*)malloc( NY * sizeof(double) );
    double* tmp = (double*)malloc( NX * sizeof(double) );

    for (int i = 0; i < (NX * NY); i++) {
        A[i] = i % 256;
    }

    for (int i = 0; i < NX; i++) {
        x[i] = i % 256;
    }

    atax_cpu_trace(A, x, y, tmp, NX, NY);
    dumpRtTmp();

    return 0;
}

