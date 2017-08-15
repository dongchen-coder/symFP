#include "./utility/rt.h"

#define M 1024
#define N 1024
#define FLOAT_N 3214212.01
#define EPS 0.005

bool varify(double* data, double* mean, double* stddev, double* symmat) {

    double* temp_mean = (double*)malloc( (M)*sizeof(double));
    double* temp_stddev = (double*)malloc( (M)*sizeof(double));
    double* temp_symmat = (double*)malloc( (M*M)*sizeof(double));

    int i, j, j1, j2;
    // Determine mean of column vectors of input data matrix
    for (j = 0; j < M; j++)
    {
        temp_mean[j] = 0.0;
        
        for (i = 0; i < N; i++)
        {
            temp_mean[j] += data[i * M + j];
        }
        
        temp_mean[j] /= (double)FLOAT_N;

        if (temp_mean[j] != mean[j]) {
            return false;
        }
    }
    
    // Determine standard deviations of column vectors of data matrix.
    for (j = 0; j < M; j++)
    {
        temp_stddev[j] = 0.0;
        
        for (i = 0; i < N; i++)
        {
            temp_stddev[j] += (data[i * M + j] - mean[j]) * (data[i * M + j] - mean[j]);
        }
        
        temp_stddev[j] /= FLOAT_N;
        temp_stddev[j] = sqrt_of_array_cell(temp_stddev, j);
        temp_stddev[j] = temp_stddev[j] <= EPS ? 1.0 : temp_stddev[j];

        if (temp_stddev[j] != stddev[j]) {
            return false;
        }
    }
    
    // Center and reduce the column vectors.
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            data[i * M + j] -= mean[j];
            data[i * M + j] /= (sqrt(FLOAT_N)*stddev[j]) ;
        }
    }
    
    // Calculate the m * m correlation matrix.
    for (j1 = 0; j1 < M-1; j1++)
    {
        temp_symmat[j1 * M + j1] = 1.0;
        
        for (j2 = j1+1; j2 < M; j2++)
        {
            temp_symmat[j1 * M + j2] = 0.0;
            
            for (i = 0; i < N; i++)
            {
                temp_symmat[j1 * M + j2] += (data[i * M + j1] * data[i * M + j2]);
            }
            
            temp_symmat[j2 * M + j1] = temp_symmat[j1 * M + j2];

            if (temp_symmat[j2 * M + j1] != symmat[j2 * M + j1]) {
                return false;
            }
        }
    }
    
    temp_symmat[(M-1) * M + M-1] = 1.0;

    return true;
}

// void correlation(int m, int n, double* data, double* mean, double* stddev, double* symmat)
// {
//     int i, j, j1, j2;
    
//     // Determine mean of column vectors of input data matrix
//     for (j = 0; j < M; j++)
//     {
//         mean[j] = 0.0;
        
//         for (i = 0; i < N; i++)
//         {
//             mean[j] += data[i * M + j];
//         }
        
//         mean[j] /= (double)FLOAT_N;
//     }
    
//     // Determine standard deviations of column vectors of data matrix.
//     for (j = 0; j < M; j++)
//     {
//         stddev[j] = 0.0;
        
//         for (i = 0; i < N; i++)
//         {
//             stddev[j] += (data[i * M + j] - mean[j]) * (data[i * M + j] - mean[j]);
//         }
        
//         stddev[j] /= FLOAT_N;
//         stddev[j] = sqrt_of_array_cell(stddev, j);
//         stddev[j] = stddev[j] <= EPS ? 1.0 : stddev[j];
//     }
    
//     // Center and reduce the column vectors.
//     for (i = 0; i < N; i++)
//     {
//         for (j = 0; j < M; j++)
//         {
//             data[i * M + j] -= mean[j];
//             data[i * M + j] /= (sqrt(FLOAT_N)*stddev[j]) ;
//         }
//     }
    
//     // Calculate the m * m correlation matrix.
//     for (j1 = 0; j1 < M-1; j1++)
//     {
//         symmat[j1 * M + j1] = 1.0;
        
//         for (j2 = j1+1; j2 < M; j2++)
//         {
//             symmat[j1 * M + j2] = 0.0;
            
//             for (i = 0; i < N; i++)
//             {
//                 symmat[j1 * M + j2] += (data[i * M + j1] * data[i * M + j2]);
//             }
            
//             symmat[j2 * M + j1] = symmat[j1 * M + j2];
//         }
//     }
    
//     symmat[(M-1) * M + M-1] = 1.0;
    
//     return;
// }

void correlation_trace(double* data, double* mean, double* stddev, double* symmat, unsigned int m, unsigned int n) {

    int i, j, j1, j2;
    
    // Determine mean of column vectors of input data matrix
    for (j = 0; j < M; j++)
    {
        mean[j] = 0.0;
        
        for (i = 0; i < N; i++)
        {
            mean[j] += data[i * M + j];
            rtTmpAccess(j + (M * N));   // load mean[j]
            rtTmpAccess(i * M + j);     // load data[i * M + j]
            rtTmpAccess(j + (M * N));   // store mean[j]
        }
        
        mean[j] /= (double)FLOAT_N;
        rtTmpAccess(j + (M * N));   // load mean[j]
        rtTmpAccess(j + (M * N));   // store mean[j]
    }
    
    // Determine standard deviations of column vectors of data matrix.
    for (j = 0; j < M; j++)
    {
        stddev[j] = 0.0;
        
        for (i = 0; i < N; i++)
        {
            stddev[j] += (data[i * M + j] - mean[j]) * (data[i * M + j] - mean[j]);
            rtTmpAccess(i * M + j);                 // load data[i * M + j];
            rtTmpAccess(j + (M * N));               // load mean[j]
            rtTmpAccess(i * M + j);                 // load data[i * M + j];
            rtTmpAccess(j + (M * N));               // load mean[j]
            rtTmpAccess(j + (M * N) + M);           // load stddev[j]
            rtTmpAccess(j + (M * N) + M);           // store stddev[j]
        }
        
        stddev[j] /= FLOAT_N;
        stddev[j] = sqrt_of_array_cell(stddev, j);
        stddev[j] = stddev[j] <= EPS ? 1.0 : stddev[j];
    }
    
    // Center and reduce the column vectors.
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            data[i * M + j] -= mean[j];
            data[i * M + j] /= (sqrt(FLOAT_N)*stddev[j]) ;
            rtTmpAccess(i * M + j);         // load data[i * M + j]
            rtTmpAccess(j + (M * N));       // load mean[j]         
            rtTmpAccess(i * M + j);         // store data[i * M + j]
            rtTmpAccess(j + (M * N + M))    // load stddev[j]
            rtTmpAccess(i * M + j);         // load data[i * M + j]
            rtTmpAccess(i * M + j)          // store data[i * M + j]
        }
    }
    
    // Calculate the m * m correlation matrix.
    for (j1 = 0; j1 < M-1; j1++)
    {
        symmat[j1 * M + j1] = 1.0;
        
        for (j2 = j1+1; j2 < M; j2++)
        {
            symmat[j1 * M + j2] = 0.0;
            
            for (i = 0; i < N; i++)
            {
                symmat[j1 * M + j2] += (data[i * M + j1] * data[i * M + j2]);
                rtTmpAccess(i * M + j1);                        // load data[i * M + j1]
                rtTmpAccess(i * M + j2);                        // load data[i * M + j2]
                rtTmpAccess(j1 * M + j2 + (M * N + 2 * M));     // load symmat[j1 * M + j2]
                rtTmpAccess(j1 * M + j2 + (M * N + 2 * M));     // store symmat[j1 * M + j2]
            }
            
            symmat[j2 * M + j1] = symmat[j1 * M + j2];
            rtTmpAccess(j1 * M + j2 + (M * N + 2 * M));     // load symmat[j1 * M + j2]
            rtTmpAccess(j2 * M + j1 + (M * N + 2 * M));     // store symmat[j1 * M + j2]
        }
    }
    
    symmat[(M-1) * M + M-1] = 1.0;
    
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

    correlation_trace(data, mean, stddev, symmat, M, N);
    dumpRtTmp();

    if (varify(data, mean, stddev, symmat)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }

    return 0;
}

