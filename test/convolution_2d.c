//void conv2D(int ni, int nj, DATA_TYPE POLYBENCH_2D(A, NI, NJ, ni, nj), DATA_TYPE POLYBENCH_2D(B, NI, NJ, ni, nj))
#define NI 1024
#define NJ 1024

void conv2D(int ni, int nj, double* A, double* B)
{
    int i, j;
    double c11, c12, c13, c21, c22, c23, c31, c32, c33;
    
    c11 = +0.2;  c21 = +0.5;  c31 = -0.8;
    c12 = -0.3;  c22 = +0.6;  c32 = -0.9;
    c13 = +0.4;  c23 = +0.7;  c33 = +0.10;
    
    
    for (i = 1; i < NI - 1; ++i) // 0
    {
        for (j = 1; j < NJ - 1; ++j) // 1
        {
            B[i * NJ + j] =  c11 * A[(i - 1) * NJ + (j - 1)]
                    +  c12 * A[(i + 0) * NJ + (j - 1)]
                    +  c13 * A[(i + 1) * NJ + (j - 1)]
                    +  c21 * A[(i - 1) * NJ + (j + 0)]
                    +  c22 * A[(i + 0) * NJ + (j + 0)]
                    +  c23 * A[(i + 1) * NJ + (j + 0)]
                    +  c31 * A[(i - 1) * NJ + (j + 1)]
                    +  c32 * A[(i + 0) * NJ + (j + 1)]
                    +  c33 * A[(i + 1) * NJ + (j + 1)];
        }
    }
}