#include "../utility/rt.h"

#define M 1024
#define N 1024
#define FLOAT_N 3214212.01
#define EPS 0.005

bool varify(double* data, double* mean, double* symmat) {

    double* temp_mean = (double*)malloc( (M*N)*sizeof(double));
    double* temp_symmat = (double*)malloc( (M*M)*sizeof(double));
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

    /* Determine mean of column vectors of input data matrix */
    for (j = 0; j < M; j++)
    {
        mean[j] = 0.0;
        for (i = 0; i < N; i++)
        {
            mean[j] = mean[j] + data[i * M + j];
            rtTmpAccess(j + (M * N));
            rtTmpAccess(i * M + j);
            rtTmpAccess(j + (M * N));
        }
        mean[j] = mean[j] / FLOAT_N;
        rtTmpAccess(j + (M * N));
        rtTmpAccess(j + (M * N));
    }
    
    /* Center the column vectors. */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            data[i * M + j] = data[i * M + j] - mean[j];
            rtTmpAccess(i * M + j);
            rtTmpAccess(j + (M * N));
            rtTmpAccess(i * M + j);
        }
    }
    
    /* Calculate the m * m covariance matrix. */
    for (j1 = 0; j1 < M; j1++)
    {
        for (j2 = j1; j2 < M; j2++)
        {
            symmat[j1 * M + j2] = 0.0;
            for (i = 0; i < N; i++)
            {
                symmat[j1 * M + j2] = symmat[j1 * M + j2] + data[i * M + j1] * data[i * M + j2];
                rtTmpAccess(i * M + j1);                    // load data[i * M + j1]
                rtTmpAccess(i * M + j2);                    // load data[i * M + j2]
                rtTmpAccess(j1 * M + j2 + (M + M * N));     // load symmat[j1 * M + j2]
                rtTmpAccess(j1 * M + j2 + (M + M * N));     // store symmat[j1 * M + j2]        
            }
            symmat[j2 * M + j1] = symmat[j1 * M + j2];
            rtTmpAccess(j1 * M + j2 + (M + M*N));
            rtTmpAccess(j2 * M + j1 + (M + M*N));
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

    covariance_trace(data, mean, symmat, M, N);
    dumpRtTmp();

    if (varify(data, mean, symmat)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }

    return 0;
}
