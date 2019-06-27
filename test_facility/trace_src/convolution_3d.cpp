#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

#ifdef ORG
	#define NI 256
	#define NJ 256
	#define NK 256
#elif defined(TX)
	#define NI 256
	#define NJ 256
	#define NK 512
#elif defined(FX)
	#define NI 256
	#define NJ 512
	#define NK 512
#elif defined(EX)
	#define NI 512
    #define NJ 512
    #define NK 512
#endif


#define A_OFFSET 0
#define B_OFFSET 

void conv3D_trace(double* A, double* B)
{
    int i, j, k;
    double c11, c12, c13, c21, c22, c23, c31, c32, c33;
    
    c11 = +2;  c21 = +5;  c31 = -8;
    c12 = -3;  c22 = +6;  c32 = -9;
    c13 = +4;  c23 = +7;  c33 = +10;

    for (int i = 1; i < NI - 1; ++i) // 0
    {
        for (int j = 1; j < NJ - 1; ++j) // 1
        {
            for (int k = 1; k < NK -1; ++k) // 2
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

                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j - 1) * NK + (k - 1));        
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j - 1) * NK + (k - 1));
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j - 1) * NK + (k - 1));
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j - 1) * NK + (k - 1));
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j - 1) * NK + (k - 1));
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j - 1) * NK + (k - 1));
                rtTmpAccess(A_OFFSET + (i + 0) * NJ * NK + (j - 1) * NK + (k + 0));
                rtTmpAccess(A_OFFSET + (i + 0) * NJ * NK + (j + 0) * NK + (k + 0));
                rtTmpAccess(A_OFFSET + (i + 0) * NJ * NK + (j + 1) * NK + (k + 0));
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j - 1) * NK + (k + 1));
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j - 1) * NK + (k + 1));
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j + 0) * NK + (k + 1));
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j + 0) * NK + (k + 1));
                rtTmpAccess(A_OFFSET + (i - 1) * NJ * NK + (j + 1) * NK + (k + 1));
                rtTmpAccess(A_OFFSET + (i + 1) * NJ * NK + (j + 1) * NK + (k + 1));
                rtTmpAccess(B_OFFSET + (i * NJ * NK + j * NK + k));
            }
        }
    }

}

int main() 
{
    double* A = (double*)malloc( (NI * NJ * NK)*sizeof(double));
    double* B = (double*)malloc( (NI * NJ * NK)*sizeof(double));

    for (int i = 0; i < (NI * NJ * NK); i++) {
            A[i] = i % 256;
    }
    
#ifdef RD
    InitRD();
#endif
    
    conv3D_trace(A, B);
    
#ifdef PROFILE_RT
    dumpRtTmp();
    RTtoMR_AET();
    dumpMR();
#endif
    
#ifdef RD
    FiniRD();
#endif

    return 0;
}
