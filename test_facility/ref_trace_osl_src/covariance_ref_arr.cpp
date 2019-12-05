
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

#define DATA_OFFSET 0
#define MEAN_OFFSET M * N
#define SYMMAT_OFFSET M * N + M

#define FLOAT_N 3214212.01
#define EPS 0.005

bool varify(double* data, double* mean, double* symmat) {

    double* temp_mean = (double*)malloc( (M*N)*sizeof(double));
    double* temp_symmat = (double*)malloc( (M*M)*sizeof(double));

	int j, i, j1, j2;
    /* Determine mean of column vectors of input data matrix */
    for (j = 0; j < M; j++)
    {
        temp_mean[j] = 0.0;
        for (i = 0; i < N; i++)
        {
            temp_mean[j] += data[i * M + j];
        }
        temp_mean[j] = temp_mean[j] / FLOAT_N;

        if (temp_mean[j] != mean[j]) 
        {
            return false;
        }
    }
    
    /* Center the column vectors. */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            data[i * M + j] = data[i * M + j] - mean[j];
        }
    }
    
    /* Calculate the m * m covariance matrix. */
    for (j1 = 0; j1 < M; j1++)
    {
        for (j2 = j1; j2 < M; j2++)
        {
            temp_symmat[j1 * M + j2] = 0.0;
            for (i = 0; i < N; i++)
            {
                temp_symmat[j1 * M + j2] = temp_symmat[j1 * M + j2] + data[i * M + j1] * data[i * M + j2];
            }
            temp_symmat[j2 * M + j1] = temp_symmat[j1 * M + j2];

            if (temp_symmat[j2 * M + j1] != symmat[j2 * M + j1]) 
            {
                return false;
            }
        }
    }

    return true;
}

// void covariance(int m, int n, double* data, double* symmat, double* mean)
// {
//     int i, j, j1,j2;
    
//     /* Determine mean of column vectors of input data matrix */
//     for (j = 0; j < M; j++)
//     {
//         mean[j] = 0.0;
//         for (i = 0; i < N; i++)
//         {
//             mean[j] += data[i * M + j];
//         }
//         mean[j] /= float_n;
//     }
    
//     /* Center the column vectors. */
//     for (i = 0; i < N; i++)
//     {
//         for (j = 0; j < M; j++)
//         {
//             data[i * M + j] -= mean[j];
//         }
//     }
    
//     /* Calculate the m * m covariance matrix. */
//     for (j1 = 0; j1 < M; j1++)
//     {
//         for (j2 = j1; j2 < M; j2++)
//         {
//             symmat[j1 * M + j2] = 0.0;
//             for (i = 0; i < N; i++)
//             {
//                 symmat[j1 * M + j2] += data[i * M + j1] * data[i * M + j2];
//             }
//             symmat[j2 * M + j1] = symmat[j1 * M + j2];
//         }
//     }
// }

void covariance_trace(double* data, double* mean, double* symmat, unsigned int m, unsigned int n) {
    
    int j, i, j1, j2;
    
    /* Determine mean of column vectors of input data matrix */
    for (j = 0; j < M; j++)
    {
        mean[j] = 0.0;
        rtTmpAccess(MEAN_OFFSET + j, 0, 0);
        
        for (i = 0; i < N; i++)
        {
            mean[j] = mean[j] + data[i * M + j];
            rtTmpAccess(DATA_OFFSET + i * M + j, 1, 1);
            rtTmpAccess(MEAN_OFFSET + j, 2, 0);
            rtTmpAccess(MEAN_OFFSET + j, 3, 0);
        }
        mean[j] = mean[j] / FLOAT_N;
        rtTmpAccess(MEAN_OFFSET + j, 4, 0);
        rtTmpAccess(MEAN_OFFSET + j, 5, 0);
    }
    
    /* Center the column vectors. */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            data[i * M + j] = data[i * M + j] - mean[j];
            rtTmpAccess(MEAN_OFFSET + j, 6, 0);
            rtTmpAccess(DATA_OFFSET + i * M + j, 7, 1);
            rtTmpAccess(DATA_OFFSET + i * M + j, 8, 1);
        }
    }
    
    /* Calculate the m * m covariance matrix. */
    for (j1 = 0; j1 < M; j1++)
    {
        for (j2 = j1; j2 < M; j2++)
        {
            symmat[j1 * M + j2] = 0.0;
            rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 9, 2);
            
            for (i = 0; i < N; i++)
            {
                symmat[j1 * M + j2] = symmat[j1 * M + j2] + data[i * M + j1] * data[i * M + j2];
                rtTmpAccess(DATA_OFFSET + i * M + j1, 10, 1);
                rtTmpAccess(DATA_OFFSET + i * M + j2, 11, 1);
                rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 12, 2);
                rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 13, 2);
            }
            symmat[j2 * M + j1] = symmat[j1 * M + j2];
            rtTmpAccess(SYMMAT_OFFSET + j1 * M + j2, 14, 2);
            rtTmpAccess(SYMMAT_OFFSET + j2 * M + j1, 15, 2);
        }
    }
    return;
}
 


int main() {
    double* data = (double*)malloc( (M*N)*sizeof(double));
    double* mean = (double*)malloc( M*sizeof(double));
    double* symmat = (double*)malloc( (M*M)*sizeof(double));

    for (int i = 0; i < M*N; i++) {
        data[i] = i % 256;
    }

#ifdef PAPI_TIMER
    PAPI_timer_init();
    PAPI_timer_start();
#endif
    covariance_trace(data, mean, symmat, M, N);
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
    /*
    if (varify(data, mean, symmat)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }
    */
    return 0;
}
