
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define NI 1024
	#define NJ 1024
#elif defined(TX)
	#define NI 1448
	#define NJ 1448
#elif defined(FX)
	#define NI 2048
	#define NJ 2048
#elif defined(EX)
	#define NI 2896
	#define NJ 2896
#endif

#define A_OFFSET 0
#define B_OFFSET NI * NJ

void conv2D_trace(double* A, double* B)
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
            
            rtTmpAccess(A_OFFSET + (i - 1) * NJ + (j - 1), 0, 0);
            rtTmpAccess(A_OFFSET + (i + 0) * NJ + (j - 1), 1, 0);
            rtTmpAccess(A_OFFSET + (i + 1) * NJ + (j - 1), 2, 0);
            rtTmpAccess(A_OFFSET + (i - 1) * NJ + (j + 0), 3, 0);
            rtTmpAccess(A_OFFSET + (i + 0) * NJ + (j + 0), 4, 0);
            rtTmpAccess(A_OFFSET + (i + 1) * NJ + (j + 0), 5, 0);
            rtTmpAccess(A_OFFSET + (i - 1) * NJ + (j + 1), 6, 0);
            rtTmpAccess(A_OFFSET + (i + 0) * NJ + (j + 1), 7, 0);
            rtTmpAccess(A_OFFSET + (i + 1) * NJ + (j + 1), 8, 0);
            rtTmpAccess(B_OFFSET + i * NJ + j, 9, 1);
            
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

    conv2D_trace(A, B);
    
    OSL_ref(0);

    return 0;
}
