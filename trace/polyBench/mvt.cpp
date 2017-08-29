#include "../utility/rt.h"

#define N 1024

bool varify(double* a, double* x1, double* x2, double* y1, double* y2)
{
    int i, j;

    double* tempX1 = (double*)malloc( N * sizeof(double));
    double* tempX2 = (double*)malloc( N * sizeof(double));
    
    for (i=0; i< N; i++)
    {
        for (j=0; j < N; j++)
        {
            //x1[i] = x1[i] + a[i][j] * y1[j];
            tempX1[i] = tempX1[i] + a[i * N + j] * y1[j];
        }
        if (tempX1[i] != x1[i])
        {
            return false;
        }
    }
    
    for (i=0; i < N; i++)
    {
        for (j=0; j< N; j++)
        {
            //x2[i] = x2[i] + a[j][i] * y2[j];
            tempX2[i] = tempX2[i] + a[j * N + i] * y2[j];
        }
        if (tempX2[i] != x2[i])
        {
            return false;
        }
    }
    return true;
}

// void runMvt(int n, double* a, double* x1, double* x2, double* y1, double* y2)
// {
//     int i, j;
    
//     for (i=0; i< N; i++)
//     {
//         for (j=0; j < N; j++)
//         {
//             //x1[i] = x1[i] + a[i][j] * y1[j];
//             x1[i] = x1[i] + a[i * N + j] * y1[j];
//         }
//     }
    
//     for (i=0; i < N; i++)
//     {
//         for (j=0; j< N; j++)
//         {
//             //x2[i] = x2[i] + a[j][i] * y2[j];
//             x2[i] = x2[i] + a[j * N + i] * y2[j];
//         }
//     }
// }

void runMvt_trace(int n, double* a, double* x1, double* x2, double* y1, double* y2)
{
    int i, j;
    
    for (i=0; i< N; i++)
    {
        for (j=0; j < N; j++)
        {
            //x1[i] = x1[i] + a[i][j] * y1[j];
            x1[i] = x1[i] + a[i * N + j] * y1[j];
            rtTmpAccess(i * N + j);     // load a
            rtTmpAccess(j + N*N);       // load y1
            rtTmpAccess(i + (N+1)*N);   // load x1
            rtTmpAccess(i + (N+1)*N);   // store x1
        }
    }
    
    for (i=0; i < N; i++)
    {
        for (j=0; j< N; j++)
        {
            //x2[i] = x2[i] + a[j][i] * y2[j];
            x2[i] = x2[i] + a[j * N + i] * y2[j];
            rtTmpAccess(j * N + i);     // load a
            rtTmpAccess(j + (N+2)*N);   // load y2
            rtTmpAccess(i + (N+3)*N);   // load x2
            rtTmpAccess(i + (N+3)*N);   // store x2
        }
    }
    return;
}

int main(int argc, char const *argv[])
{
    double* a = (double*)malloc( (N*N) * sizeof(double));
    double* y1 = (double*)malloc( N * sizeof(double));
    double* x1 = (double*)malloc( N * sizeof(double));
    double* y2 = (double*)malloc( N * sizeof(double));
    double* x2 = (double*)malloc( N * sizeof(double));
    

    for (int i = 0; i < N; ++i)
    {
        y1[i] = i % 256;
        y2[i] = i / 100;
    }

    for (int i = 0; i < N*N; ++i)
    {
        a[i] = i / 10;
    }

    runMvt_trace(N, a, x1, x2, y1, y2);
    dumpRtTmp();

    if (varify(alpha, beta, a, x1, x2, y1, y2)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }
    return 0;
}
