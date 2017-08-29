#include "../utility/rt.h"
#define N 1024

bool varify(double alpha, double beta, double* A, double* B, double* tmp, double* x, double* y)
{
    int i, j;

    double* tempY = (double*)malloc( N * sizeof(double));
    
    for (i = 0; i < N; i++)
    {
        tmp[i] = 0;
        tempY[i] = 0;
        for (j = 0; j < N; j++)
        {
            tmp[i] = A[i * N + j] * x[j] + tmp[i];
            tempY[i] = B[i * N + j] * x[j] + tempY[i];
        }
        
        tempY[i] = alpha * tmp[i] + beta * tempY[i];
        if (tempY[i] != y[i])
        {
            return false;
        }
    }
    
    free(tempY);
    return true;
}

// void gesummv(int n, double alpha, double beta, double* A, double* B, double* tmp, double* x, double* y)
// {
//     int i, j;
    
//     for (i = 0; i < N; i++)
//     {
//         tmp[i] = 0;
//         y[i] = 0;
//         for (j = 0; j < N; j++)
//         {
//             tmp[i] = A[i * N + j] * x[j] + tmp[i];
//             y[i] = B[i * N + j] * x[j] + y[i];
//         }
        
//         y[i] = alpha * tmp[i] + beta * y[i];
//     }
// }

void gesummv_trace(int n, double alpha, double beta, double* A, double* B, double* tmp, double* x, double* y)
{
    int i, j;
    
    for (i = 0; i < N; i++)
    {
        tmp[i] = 0;
        y[i] = 0;
        for (j = 0; j < N; j++)
        {
            tmp[i] = A[i * N + j] * x[j] + tmp[i];
            rtTmpAccess(i * N + j);         // load A
            rtTmpAccess(j + N*N);           // load x
            rtTmpAccess(i + (N+1)*N );      // load tmp
            rtTmpAccess(i + (N+1)*N );      // store tmp
            
            y[i] = B[i * N + j] * x[j] + y[i];
            rtTmpAccess(i * N + j + (N+2)*N);  // load B
            rtTmpAccess(j + N*N);              // load x
            rtTmpAccess(i + (N*2 + 2)*N);      // load y
            rtTmpAccess(i + (N*2 + 2)*N);      // store y
        }
        
        y[i] = alpha * tmp[i] + beta * y[i];
        rtTmpAccess(i + (N+1)*N);       // load tmp
        rtTmpAccess(i + (N*2 + 2)*N);   // load y
        rtTmpAccess(i + (N*2 + 2)*N);   // store y
    }
    return;
}

int main(int argc, char const *argv[])
{
    double* A = (double*)malloc( (N*N) * sizeof(double));
    double* x = (double*)malloc( N * sizeof(double));
    double* tmp = (double*)malloc( N * sizeof(double));
    double* B = (double*)malloc( (N*N) * sizeof(double));
    double* y = (double*)malloc( N * sizeof(double));

    for (int i = 0; i < N; ++i)
    {
        x[i] = i % 256;
    }

    for (int i = 0; i < N*N; ++i)
    {
        A[i] = i / 10;
        B[i] = i / 25;
    }

    double alpha = 1.0;
    double beta = 1.5;

    gesummv_trace(N, alpha, beta, A, B, tmp, x, y);
    dumpRtTmp();

    if (varify(alpha, beta, A, B, tmp, x, y)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }
    return 0;
}
