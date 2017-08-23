#include "utility/rt.h"

#define NI 1024
#define NJ 1024

bool varify(double alpha, double beta, double* A, double* C)
{
    int i, j, k;
    
    /*  C := alpha*A*A' + beta*C */
    double tempC = 1.0;
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            tempC *= beta;
        }
    }
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            for (k = 0; k < NJ; k++)
            {
                tempC += alpha * A[i * NJ + k] * A[j * NJ + k];
            }
            if (tempC != C[i * NI + j]) {
                return false;
            }
        }
    }

    return true;
}

// void syrk(int ni, int nj, double alpha, double beta, double* A, double* C)
// {
//     int i, j, k;
    
//     /*  C := alpha*A*A' + beta*C */
//     for (i = 0; i < NI; i++)
//     {
//         for (j = 0; j < NI; j++)
//         {
//             C[i * NI + j] *= beta;
//         }
//     }
    
//     for (i = 0; i < NI; i++)
//     {
//         for (j = 0; j < NI; j++)
//         {
//             for (k = 0; k < NJ; k++)
//             {
//                 C[i * NI + j] += alpha * A[i * NJ + k] * A[j * NJ + k];
//             }
//         }
//     }
// }

void syrk_trace(int ni, int nj, double alpha, double beta, double* A, double* C)
{
    int i, j, k;
    
    /*  C := alpha*A*A' + beta*C */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            C[i * NI + j] = C[i * NI + j] * beta;
            rtTmpAccess(i * NI + j);        // load C
            rtTmpAccess(i * NI + j);        // store C
        }
    }
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NI; j++)
        {
            for (k = 0; k < NJ; k++)
            {
                C[i * NI + j] = C[i * NI + j] + alpha * A[i * NJ + k] * A[j * NJ + k];
                rtTmpAccess(i * NJ + k + (N * N));      // load A
                rtTmpAccess(j * NJ + k + (N * N));      // load A
                rtTmpAccess(i * NI + j);                // load C
                rtTmpAccess(i * NI + j);                // store C
            }
        }
    }
    return;
}

int main(int argc, char const *argv[])
{
    double* A = (double*)malloc( (N*N) * sizeof(double));
    double* C = (double*)malloc( (N*N) * sizeof(double));

    for (int i = 0; i < N*N; ++i)
    {
        A[i] = i / 10;
        C[i] = 1.0;
    }

    double alpah = 0.0;
    double beta = 1.5;

    syrk_trace(NI, NJ, alpha, beta, A, C);
    dumpRtTmp();

    if (varify(alpha, beta, A, C)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }
    return 0;
}