#include "../utility/rt.h"

#define N 1024

bool varify(double alpha, double beta, double* A, double* u1, double* v1, double* u2, double* v2, double* w, double* x, double* y, double* z)
{
    int i,j;

    double* tempX = (double*)malloc( N * sizeof(double));
    double* tempW = (double*)malloc( N * sizeof(double));
    
    for (i = 0; i < N; i++)
    {
        double tempA = 0.0;
        for (j = 0; j < N; j++)
        {
            tempA[i * N + j] = tempA + u1[i] * v1[j] + u2[i] * v2[j];
            if (A[i * N + j] != tempA)
            {
                return false;
            }
        }
    }
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            tempX[i] = tempX[i] + beta * A[j * N + i] * y[j];
        }
        if (tempX[i] != x[i])
        {
            return false;
        }
    }
    
    for (i = 0; i < N; i++)
    {
        x[i] = x[i] + z[i];
    }
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            tempW[i] = tempW[i] +  alpha * A[i * N + j] * x[j];
        }
        if (tempW[i] != w[i]) 
        {
            return false;
        }
    }

    free(tempW);
    free(tempX);

    return true;

}

// void gemver(int n, double alpha, double beta, double* A, double* u1, double* v1, double* u2, double* v2, double* w, double* x, double* y, double* z)
// {
//     int i,j;
    
//     for (i = 0; i < N; i++)
//     {
//         for (j = 0; j < N; j++)
//         {
//             A[i * N + j] = A[i * N + j] + u1[i] * v1[j] + u2[i] * v2[j];
//         }
//     }
    
//     for (i = 0; i < N; i++)
//     {
//         for (j = 0; j < N; j++)
//         {
//             x[i] = x[i] + beta * A[j * N + i] * y[j];
//         }
//     }
    
//     for (i = 0; i < N; i++)
//     {
//         x[i] = x[i] + z[i];
//     }
    
//     for (i = 0; i < N; i++)
//     {
//         for (j = 0; j < N; j++)
//         {
//             w[i] = w[i] +  alpha * A[i * N + j] * x[j];
//         }
//     }
// }

void gemver_trace(int n, double alpha, double beta, double* A, double* u1, double* v1, double* u2, double* v2, double* w, double* x, double* y, double* z)
{
    int i,j;
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i * N + j] = A[i * N + j] + u1[i] * v1[j] + u2[i] * v2[j];
            rtTmpAccess(i);                 // load u1
            rtTmpAccess(j + N);             // load v1
            rtTmpAccess(i + 2*N);           // load u2
            rtTmpAccess(j + 3*N);           // load v2
            rtTmpAccess(i * N + j + 4*N);   // load A
            rtTmpAccess(i * N + j + 4*N);   // store A
        }
    }
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            x[i] = x[i] + beta * A[j * N + i] * y[j];
            rtTmpAccess(j * N + i + 4*N);   // load A
            rtTmpAccess(j + (4+N)*N);       // load y
            rtTmpAccess(i + (5+N)*N);       // load x
            rtTmpAccess(i + (5+N)*N);       // store x
        }
    }
    
    for (i = 0; i < N; i++)
    {
        x[i] = x[i] + z[i];
        rtTmpAccess(i + (5+N)*N);       // load x
        rtTmpAccess(i + (6*N)*N);       // load z
        rtTmpAccess(i + (5+N)*N);       // store x
    }
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            w[i] = w[i] +  alpha * A[i * N + j] * x[j];
            rtTmpAccess(i * N + j + 4*N);   // load A
            rtTmpAccess(j + (5+N)*N);       // load x
            rtTmpAccess(i + (6*N)*N);       // load w
            rtTmpAccess(i + (6*N)*N);       // store w

        }
    }
    return;
}

int main(int argc, char const *argv[])
{
    double* u1 = (double*)malloc( N * sizeof(double));
    double* v1 = (double*)malloc( N * sizeof(double));
    double* u2 = (double*)malloc( N * sizeof(double));
    double* v2 = (double*)malloc( N * sizeof(double));
    double* A = (double*)malloc( (N*N) * sizeof(double));
    double* y = (double*)malloc( N * sizeof(double));
    double* x = (double*)malloc( N * sizeof(double));
    double* z = (double*)malloc( N * sizeof(double));
    double* w = (double*)malloc( N * sizeof(double));

    for (int i = 0; i < N; ++i)
    {
        u1[i] = i % 256;
        u2[i] = i % 32;
        v1[i] = i % 2;
        v2[i] = i % 512;
        y[i] = i % 4;
        z[i] = i % 64;
        x[i] = 0.0;
        w[i] = 0.0;
    }

    for (int i = 0; i < N*N; ++i)
    {
        A[i] = 0.0;
    }

    double alpha = 1.0;
    double beta = 1.5;

    gemver_trace(N, alpha, beta, A, u1, v1, u2, v2, w, x, y, z);
    dumpRtTmp();

    if (varify(alpha, beta, A, u1, v1, u2, v2, w, x, y, z)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }
    return 0;
}
