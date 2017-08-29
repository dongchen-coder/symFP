#include "../utility/rt.h"

#define NI 1024
#define NJ 1024

bool varify(double* A, double* B) 
{
    double c11, c12, c13, c21, c22, c23, c31, c32, c33;
    
    c11 = +2;  c21 = +5;  c31 = -8;
    c12 = -3;  c22 = +6;  c32 = -9;
    c13 = +4;  c23 = +7;  c33 = +10;

    for (i = 1; i < NI - 1; ++i) // 0
    {
        for (j = 1; j < NJ - 1; ++j) // 1
        {
            for (k = 1; k < NK -1; ++k) // 2
            {
                /*
                B[i][j][k] = c11 * A[(i - 1)][(j - 1)][(k - 1)]  +  c13 * A[(i + 1)][(j - 1)][(k - 1)]
                +   c21 * A[(i - 1)][(j - 1)][(k - 1)]  +  c23 * A[(i + 1)][(j - 1)][(k - 1)]
                +   c31 * A[(i - 1)][(j - 1)][(k - 1)]  +  c33 * A[(i + 1)][(j - 1)][(k - 1)]
                +   c12 * A[(i + 0)][(j - 1)][(k + 0)]  +  c22 * A[(i + 0)][(j + 0)][(k + 0)]
                +   c32 * A[(i + 0)][(j + 1)][(k + 0)]  +  c11 * A[(i - 1)][(j - 1)][(k + 1)]
                +   c13 * A[(i + 1)][(j - 1)][(k + 1)]  +  c21 * A[(i - 1)][(j + 0)][(k + 1)]
                +   c23 * A[(i + 1)][(j + 0)][(k + 1)]  +  c31 * A[(i - 1)][(j + 1)][(k + 1)]
                +   c33 * A[(i + 1)][(j + 1)][(k + 1)];
                 */
                double temp = c11 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                            +   c13 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                            +   c21 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                            +   c23 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                            +   c31 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                            +   c33 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k - 1)]
                            +   c12 * A[(i + 0) * NJ * NK + (j - 1) * NK + (k + 0)]
                            +   c22 * A[(i + 0) * NJ * NK + (j + 0) * NK + (k + 0)]
                            +   c32 * A[(i + 0) * NJ * NK + (j + 1) * NK + (k + 0)]
                            +   c11 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k + 1)]
                            +   c13 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k + 1)]
                            +   c21 * A[(i - 1) * NJ * NK + (j + 0) * NK + (k + 1)]
                            +   c23 * A[(i + 1) * NJ * NK + (j + 0) * NK + (k + 1)]
                            +   c31 * A[(i - 1) * NJ * NK + (j + 1) * NK + (k + 1)]
                            +   c33 * A[(i + 1) * NJ * NK + (j + 1) * NK + (k + 1)];
                if (temp != B[i * NJ * NK + j * NK + k]) 
                    return false;
            }
        }
    }

    return true;

}

// void conv3D(int ni, int nj, int nk, double* A, double* B)
// {
//     int i, j, k;
//     double c11, c12, c13, c21, c22, c23, c31, c32, c33;
    
//     c11 = +2;  c21 = +5;  c31 = -8;
//     c12 = -3;  c22 = +6;  c32 = -9;
//     c13 = +4;  c23 = +7;  c33 = +10;
    
//     for (i = 1; i < NI - 1; ++i) // 0
//     {
//         for (j = 1; j < NJ - 1; ++j) // 1
//         {
//             for (k = 1; k < NK -1; ++k) // 2
//             {
//                 /*
//                 B[i][j][k] = c11 * A[(i - 1)][(j - 1)][(k - 1)]  +  c13 * A[(i + 1)][(j - 1)][(k - 1)]
//                 +   c21 * A[(i - 1)][(j - 1)][(k - 1)]  +  c23 * A[(i + 1)][(j - 1)][(k - 1)]
//                 +   c31 * A[(i - 1)][(j - 1)][(k - 1)]  +  c33 * A[(i + 1)][(j - 1)][(k - 1)]
//                 +   c12 * A[(i + 0)][(j - 1)][(k + 0)]  +  c22 * A[(i + 0)][(j + 0)][(k + 0)]
//                 +   c32 * A[(i + 0)][(j + 1)][(k + 0)]  +  c11 * A[(i - 1)][(j - 1)][(k + 1)]
//                 +   c13 * A[(i + 1)][(j - 1)][(k + 1)]  +  c21 * A[(i - 1)][(j + 0)][(k + 1)]
//                 +   c23 * A[(i + 1)][(j + 0)][(k + 1)]  +  c31 * A[(i - 1)][(j + 1)][(k + 1)]
//                 +   c33 * A[(i + 1)][(j + 1)][(k + 1)];
//                  */
//                 B[i * NJ * NK + j * NK + k] = c11 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
//                                             +   c13 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k - 1)]
//                                             +   c21 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
//                                             +   c23 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k - 1)]
//                                             +   c31 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
//                                             +   c33 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k - 1)]
//                                             +   c12 * A[(i + 0) * NJ * NK + (j - 1) * NK + (k + 0)]
//                                             +   c22 * A[(i + 0) * NJ * NK + (j + 0) * NK + (k + 0)]
//                                             +   c32 * A[(i + 0) * NJ * NK + (j + 1) * NK + (k + 0)]
//                                             +   c11 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k + 1)]
//                                             +   c13 * A[(i + 1) * NJ * NK + (j - 1) * NK + (k + 1)]
//                                             +   c21 * A[(i - 1) * NJ * NK + (j + 0) * NK + (k + 1)]
//                                             +   c23 * A[(i + 1) * NJ * NK + (j + 0) * NK + (k + 1)]
//                                             +   c31 * A[(i - 1) * NJ * NK + (j + 1) * NK + (k + 1)]
//                                             +   c33 * A[(i + 1) * NJ * NK + (j + 1) * NK + (k + 1)];
//             }
//         }
//     }
// }


void conv3D_trace(double* A, double* B, int ni, int nj, int nk)
{
    int i, j, k;
    double c11, c12, c13, c21, c22, c23, c31, c32, c33;
    
    c11 = +2;  c21 = +5;  c31 = -8;
    c12 = -3;  c22 = +6;  c32 = -9;
    c13 = +4;  c23 = +7;  c33 = +10;
    
    
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

        }
    }

    for (i = 1; i < NI - 1; ++i) // 0
    {
        for (j = 1; j < NJ - 1; ++j) // 1
        {
            for (k = 1; k < NK -1; ++k) // 2
            {
                /*
                B[i][j][k] = c11 * A[(i - 1)][(j - 1)][(k - 1)]  +  c13 * A[(i + 1)][(j - 1)][(k - 1)]
                +   c21 * A[(i - 1)][(j - 1)][(k - 1)]  +  c23 * A[(i + 1)][(j - 1)][(k - 1)]
                +   c31 * A[(i - 1)][(j - 1)][(k - 1)]  +  c33 * A[(i + 1)][(j - 1)][(k - 1)]
                +   c12 * A[(i + 0)][(j - 1)][(k + 0)]  +  c22 * A[(i + 0)][(j + 0)][(k + 0)]
                +   c32 * A[(i + 0)][(j + 1)][(k + 0)]  +  c11 * A[(i - 1)][(j - 1)][(k + 1)]
                +   c13 * A[(i + 1)][(j - 1)][(k + 1)]  +  c21 * A[(i - 1)][(j + 0)][(k + 1)]
                +   c23 * A[(i + 1)][(j + 0)][(k + 1)]  +  c31 * A[(i - 1)][(j + 1)][(k + 1)]
                +   c33 * A[(i + 1)][(j + 1)][(k + 1)];
                 */
                B[i * nj * nk + j * nk + k] = c11 * A[(i - 1) * nj * nk + (j - 1) * nk + (k - 1)]
                                            +   c13 * A[(i + 1) * nj * nk + (j - 1) * nk + (k - 1)]
                                            +   c21 * A[(i - 1) * nj * nk + (j - 1) * nk + (k - 1)]
                                            +   c23 * A[(i + 1) * nj * nk + (j - 1) * nk + (k - 1)]
                                            +   c31 * A[(i - 1) * nj * nk + (j - 1) * nk + (k - 1)]
                                            +   c33 * A[(i + 1) * nj * nk + (j - 1) * nk + (k - 1)]
                                            +   c12 * A[(i + 0) * nj * nk + (j - 1) * nk + (k + 0)]
                                            +   c22 * A[(i + 0) * nj * nk + (j + 0) * nk + (k + 0)]
                                            +   c32 * A[(i + 0) * nj * nk + (j + 1) * nk + (k + 0)]
                                            +   c11 * A[(i - 1) * nj * nk + (j - 1) * nk + (k + 1)]
                                            +   c13 * A[(i + 1) * nj * nk + (j - 1) * nk + (k + 1)]
                                            +   c21 * A[(i - 1) * nj * nk + (j + 0) * nk + (k + 1)]
                                            +   c23 * A[(i + 1) * nj * nk + (j + 0) * nk + (k + 1)]
                                            +   c31 * A[(i - 1) * nj * nk + (j + 1) * nk + (k + 1)]
                                            +   c33 * A[(i + 1) * nj * nk + (j + 1) * nk + (k + 1)];

                rtTmpAccess((i - 1) * nj * nk + (j - 1) * nk + (k - 1));        
                rtTmpAccess((i + 1) * nj * nk + (j - 1) * nk + (k - 1));
                rtTmpAccess((i - 1) * nj * nk + (j - 1) * nk + (k - 1));
                rtTmpAccess((i + 1) * nj * nk + (j - 1) * nk + (k - 1));
                rtTmpAccess((i - 1) * nj * nk + (j - 1) * nk + (k - 1));
                rtTmpAccess((i + 1) * nj * nk + (j - 1) * nk + (k - 1));
                rtTmpAccess((i + 0) * nj * nk + (j - 1) * nk + (k + 0));
                rtTmpAccess((i + 0) * nj * nk + (j + 0) * nk + (k + 0));
                rtTmpAccess((i + 0) * nj * nk + (j + 1) * nk + (k + 0));
                rtTmpAccess((i - 1) * nj * nk + (j - 1) * nk + (k + 1));
                rtTmpAccess((i + 1) * nj * nk + (j - 1) * nk + (k + 1));
                rtTmpAccess((i - 1) * nj * nk + (j + 0) * nk + (k + 1));
                rtTmpAccess((i + 1) * nj * nk + (j + 0) * nk + (k + 1));
                rtTmpAccess((i - 1) * nj * nk + (j + 1) * nk + (k + 1));
                rtTmpAccess((i + 1) * nj * nk + (j + 1) * nk + (k + 1));
                rtTmpAccess((i * nj * nk + j * nk + k) + (ni * nj * nk));
            }
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
