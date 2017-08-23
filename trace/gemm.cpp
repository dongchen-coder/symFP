#include "./utility/rt.h"

#define NI 1024
#define NJ 1024
#define NK 1024

void varify(double alpha, double beta, double* A, double* B, double* C) {
    int i,j,k;

    double *tempC = (double*)malloc( (NI * NJ) * sizeof(double));
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            tempC[i * NJ + j] *= beta;
            
            for (k = 0; k < NK; k++)
            {
                tempC[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
            }
            if (tempC[i * NJ + j] != C[i * NJ + j]) {
                return false;
            }
        }
    } 
    return true;  
}

// void gemm(int ni, int nj, int nk, double alpha, double beta, double* A, double* B, double* C)
// {
//     int i,j,k;
    
//     for (i = 0; i < NI; i++)
//     {
//         for (j = 0; j < NJ; j++)
//         {
//             C[i * NJ + j] *= beta;
            
//             for (k = 0; k < NK; k++)
//             {
//                 C[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
//             }
//         }
//     }
// }

void gemm_trace(int ni, int nj, int nk, double alpha, double beta, double* A, double* B, double* C) {
    int i,j,k;
    
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            C[i * NJ + j] *= beta;
            
            for (k = 0; k < NK; k++)
            {
                C[i * NJ + j] += alpha * A[i * NK + k] * B[k * NJ + j];
                rtTmpAccess(i * NK + k);                            // load A
                rtTmpAccess(k * NJ + j + (NI * NK));                // load B
                rtTmpAccess(i * NJ + j + (NI * NK) + (NK * NJ));    // load C
                rtTmpAccess(i * NJ + j + (NI * NK) + (NK * NJ));    // store C
            }
        }
    }
    return;
}

int main()
{
    double* A = (double*)malloc(NI * NK * sizeof(double));
    double* B = (double*)malloc(NK * NJ * sizeof(double));
    double* C = (double*)malloc(NI * NJ * sizeof(double));

    for (int i = 0; i < NI*NK; ++i)
    {
        A[i] = i % 256;
    }

    for (int i = 0; i < NK * NJ; ++i)
    {
        B[i] = i % 48;
    }

    double alpha = 1.0;
    double beta = 1.5;

    gemm_trace(NI, NJ, NK, alpha, beta, A, B, C);
    dumpRtTmp();

    if (varify(alpha, beta, A, B, C)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }
    return 0;
}