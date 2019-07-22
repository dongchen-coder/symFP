
#include "../utility/rt.h"
#include "../utility/data_size.h"
#include<math.h>

#ifdef ORG
	#define M 1024
	#define N 1024
#elif defined(TX)
#elif defined(FX)
#endif


#define FLOAT_N 3214212.01
#define EPS 0.005

#define DATA_OFFSET 0
#define MEAN_OFFSET M * N
#define STDDEV_OFFSET M * N + M
#define SYMMAT_OFFSET M * N + M + M

void correlation_trace(double* data, double* mean, double* stddev, double* symmat) {

	int i, j, j1, j2;

    // Determine mean of column vectors of input data matrix
    for (j = 0; j < M; j++)
    {
        mean[j] = 0.0;
		rtTmpAccess(MEAN_OFFSET + j, 0, 0);

        for (i = 0; i < N; i++)
        {
            mean[j] += data[i * M + j];
            rtTmpAccess(DATA_OFFSET + i * M + j, 1, 1);
        	rtTmpAccess(MEAN_OFFSET + j, 2, 0);
			rtTmpAccess(MEAN_OFFSET + j, 3, 0);
		}

        mean[j] /= (double)FLOAT_N;
    	rtTmpAccess(MEAN_OFFSET + j, 4, 0);
		rtTmpAccess(MEAN_OFFSET + j, 5, 0);
	}

    // Determine standard deviations of column vectors of data matrix.
    for (j = 0; j < M; j++)
    {
        stddev[j] = 0.0;
		rtTmpAccess(STDDEV_OFFSET + j, 6, 2);

        for (i = 0; i < N; i++)
        {
            stddev[j] += (data[i * M + j] - mean[j]) * (data[i * M + j] - mean[j]);
        	
			rtTmpAccess(DATA_OFFSET + i * M + j, 7, 1);
			rtTmpAccess(MEAN_OFFSET + j, 8, 0);
			rtTmpAccess(DATA_OFFSET + i * M + j, 9, 1);
			rtTmpAccess(MEAN_OFFSET + j, 10, 0);
            rtTmpAccess(STDDEV_OFFSET + j, 11, 2);
			rtTmpAccess(STDDEV_OFFSET + j, 12, 2);
		}

        stddev[j] /= FLOAT_N;
//        stddev[j] = sqrt_of_array_cell(stddev, j);
		stddev[j] = sqrt(stddev[j]);

        stddev[j] = stddev[j] <= EPS ? 1.0 : stddev[j];
    	rtTmpAccess(STDDEV_OFFSET + j, 13, 2);
		rtTmpAccess(STDDEV_OFFSET + j, 14, 2);
		rtTmpAccess(STDDEV_OFFSET + j, 15, 2);
		rtTmpAccess(STDDEV_OFFSET + j, 16, 2);
        if (stddev[j] > EPS) {
            rtTmpAccess(STDDEV_OFFSET + j, 17, 2);
        }
		rtTmpAccess(STDDEV_OFFSET + j, 18, 2);
	}

	// Center and reduce the column vectors.
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            data[i * M + j] -= mean[j];
            data[i * M + j] /= (sqrt(FLOAT_N)*stddev[j]) ;
            
            rtTmpAccess(MEAN_OFFSET + j, 19, 0);
			rtTmpAccess(DATA_OFFSET + i * M + j, 20, 1);
			rtTmpAccess(DATA_OFFSET + i * M + j, 21, 1);
            rtTmpAccess(STDDEV_OFFSET + j, 22, 2);
			rtTmpAccess(DATA_OFFSET + i * M + j, 23, 1);
			rtTmpAccess(DATA_OFFSET + i * M + j, 24, 1);
        }
    }

    // Calculate the m * m correlation matrix.
    for (j1 = 0; j1 < M-1; j1++)
    {
        symmat[j1 * M + j1] = 1.0;
		rtTmpAccess(SYMMAT_OFFSET + j1 * M + j1, 25, 3);

        for (j2 = j1+1; j2 < M; j2++)
        {
            symmat[j1 * M + j2] = 0.0;
			rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 26, 3);

            for (i = 0; i < N; i++)
            {
                symmat[j1 * M + j2] += (data[i * M + j1] * data[i * M + j2]);
                rtTmpAccess(DATA_OFFSET + i * M + j1, 27, 1);
                rtTmpAccess(DATA_OFFSET + i * M + j2, 28, 1);
                rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 29, 3);
				rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 30, 3);
            }

            symmat[j2 * M + j1] = symmat[j1 * M + j2];
			rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 31, 3);
			rtTmpAccess(SYMMAT_OFFSET + j2 * M + j1, 32, 3);
        }
    }

    symmat[(M-1) * M + M-1] = 1.0;	
	rtTmpAccess(SYMMAT_OFFSET + (M-1) * M + M-1, 33, 3);

    return;
}
 


int main() {
    double* data = (double*)malloc( (M*N)*sizeof(double));
    double* mean = (double*)malloc( M*sizeof(double));
    double* stddev = (double*)malloc( M*sizeof(double));
    double* symmat = (double*)malloc( (M*M)*sizeof(double));

    for (int i = 0; i < M; i++) {
        data[i] = i % 256;
    }

    correlation_trace(data, mean, stddev, symmat);
    
    dumpRI();

    return 0;
}
