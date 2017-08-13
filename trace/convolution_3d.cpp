#include "utility/rt.h"

#define NI 1024
#define NJ 1024
#define NK 1024

bool varify(double* A, double* B) 
{
    double c11, c12, c13, c21, c22, c23, c31, c32, c33;
    
    c11 = +0.2;  c21 = +0.5;  c31 = -0.8;
    c12 = -0.3;  c22 = +0.6;  c32 = -0.9;
    c13 = +0.4;  c23 = +0.7;  c33 = +0.10;

    for (i = 1; i < NI - 1; ++i) // 0
    {
        for (j = 1; j < NJ - 1; ++j) // 1
        {
            double temp =  c11 * A[(i - 1) * NJ + (j - 1)]
                        +  c12 * A[(i + 0) * NJ + (j - 1)]
                        +  c13 * A[(i + 1) * NJ + (j - 1)]
                        +  c21 * A[(i - 1) * NJ + (j + 0)]
                        +  c22 * A[(i + 0) * NJ + (j + 0)]
                        +  c23 * A[(i + 1) * NJ + (j + 0)]
                        +  c31 * A[(i - 1) * NJ + (j + 1)]
                        +  c32 * A[(i + 0) * NJ + (j + 1)]
                        +  c33 * A[(i + 1) * NJ + (j + 1)];
            if (B[i * NJ + j] != temp) 
                return false;
        }
    }

    return true;

}

// void conv2D(int ni, int nj, double* A, double* B)
// {
//     int i, j;
//     double c11, c12, c13, c21, c22, c23, c31, c32, c33;
    
//     c11 = +0.2;  c21 = +0.5;  c31 = -0.8;
//     c12 = -0.3;  c22 = +0.6;  c32 = -0.9;
//     c13 = +0.4;  c23 = +0.7;  c33 = +0.10;
    
    
//     for (i = 1; i < NI - 1; ++i) // 0
//     {
//         for (j = 1; j < NJ - 1; ++j) // 1
//         {
//             B[i * NJ + j] =  c11 * A[(i - 1) * NJ + (j - 1)]
//                     +  c12 * A[(i + 0) * NJ + (j - 1)]
//                     +  c13 * A[(i + 1) * NJ + (j - 1)]
//                     +  c21 * A[(i - 1) * NJ + (j + 0)]
//                     +  c22 * A[(i + 0) * NJ + (j + 0)]
//                     +  c23 * A[(i + 1) * NJ + (j + 0)]
//                     +  c31 * A[(i - 1) * NJ + (j + 1)]
//                     +  c32 * A[(i + 0) * NJ + (j + 1)]
//                     +  c33 * A[(i + 1) * NJ + (j + 1)];
//         }
//     }
// }

void conv2D_trace(double* A, double* B, int ni, int nj)
{
    int i, j;
    double c11, c12, c13, c21, c22, c23, c31, c32, c33;
    
    c11 = +0.2;  c21 = +0.5;  c31 = -0.8;
    c12 = -0.3;  c22 = +0.6;  c32 = -0.9;
    c13 = +0.4;  c23 = +0.7;  c33 = +0.10;
    
    
    for (i = 1; i < ni - 1; ++i) // 0
    {
        for (j = 1; j < nj - 1; ++j) // 1
        {
            B[i * nj + j] =  c11 * A[(i - 1) * nj + (j - 1)]
                    +  c12 * A[(i + 0) * nj + (j - 1)]
                    +  c13 * A[(i + 1) * nj + (j - 1)]
                    +  c21 * A[(i - 1) * nj + (j + 0)]
                    +  c22 * A[(i + 0) * nj + (j + 0)]
                    +  c23 * A[(i + 1) * nj + (j + 0)]
                    +  c31 * A[(i - 1) * nj + (j + 1)]
                    +  c32 * A[(i + 0) * nj + (j + 1)]
                    +  c33 * A[(i + 1) * nj + (j + 1)];

            rtTmpAccess((i - 1) * nj + (j - 1));        
            rtTmpAccess((i + 0) * nj + (j - 1));
            rtTmpAccess((i + 1) * nj + (j - 1));
            rtTmpAccess((i - 1) * nj + (j + 0));
            rtTmpAccess((i + 0) * nj + (j + 0));
            rtTmpAccess((i + 1) * nj + (j + 0));
            rtTmpAccess((i - 1) * nj + (j + 1));
            rtTmpAccess((i + 0) * nj + (j + 1));
            rtTmpAccess((i + 1) * nj + (j + 1));
            rtTmpAccess(i * nj + j + (ni * nj));

        }
    }
}

int main() 
{
    double* A = (double*)malloc( (NI * NJ)*sizeof(double));
    double* B = (double*)malloc( (NI * NJ)*sizeof(double));

    for (int i = 0; i < (NI * NJ); i++) {
            A[i] = i % 256;
    }

    conv2D_trace(A, B, NI, NJ);
    dumpRtTmp();

    if (varify(A, B)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }

    return 0;
    return 0;
}