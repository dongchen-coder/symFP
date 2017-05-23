//void conv3D(int ni, int nj, int nk, DATA_TYPE POLYBENCH_3D(A, NI, NJ, NK, ni, nj, nk), DATA_TYPE POLYBENCH_3D(B, NI, NJ, NK, ni, nj, nk))

#define NI 1024
#define NJ 1024
#define NK 1024

void conv3D(int ni, int nj, int nk, double* A, double* B)
{
    int i, j, k;
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
                B[i * NJ * NK + j * NK + k] = c11 * A[(i - 1) * NJ * NK + (j - 1) * NK + (k - 1)]
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
            }
        }
    }
}
