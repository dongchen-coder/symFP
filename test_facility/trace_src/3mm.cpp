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
	#define NL 256
	#define NK 256
	#define NM 256
#elif defined (TX)
	#define NI 512
    #define NJ 256
    #define NL 512
    #define NK 256
    #define NM 256
#elif defined (FX)
	#define NI 1024
    #define NJ 256
    #define NL 1024
    #define NK 256
    #define NM 256
#elif defined (EX)
	#define NI 2048
	#define NJ 256
	#define NL 2048
	#define NK 256
	#define NM 256
#endif

#define A_OFFSET 0
#define B_OFFSET NK * NI
#define C_OFFSET NK * NI + NK * NJ
#define D_OFFSET NK * NI + NK * NJ + NJ * NM
#define E_OFFSET NK * NI + NK * NJ + NJ * NM + NM * NL
#define F_OFFSET NK * NI + NK * NJ + NJ * NM + NM * NL + NI * NJ
#define G_OFFSET NK * NI + NK * NJ + NJ * NM + NM * NL + NI * NJ + NJ * NL


void mm3_cpu_trace(int *A, int *B, int *C, int *D, int *E, int *F, int *G) {
    
	int i, j, k;

    /* E := A*B */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            E[i * NJ + j] = 0;
			
			rtTmpAccess(E_OFFSET + i * NJ + j);
	
            for (k = 0; k < NK; ++k)
            {
                E[i * NJ + j] += A[i * NK + k] * B[k * NJ + j];
				rtTmpAccess(A_OFFSET + i * NK + k);
				rtTmpAccess(B_OFFSET + k * NJ + j);
				rtTmpAccess(E_OFFSET + i * NJ + j);
				rtTmpAccess(E_OFFSET + i * NJ + j);
            }
        }
    }

    /* F := C*D */
    for (i = 0; i < NJ; i++)
    {
        for (j = 0; j < NL; j++)
        {
            F[i * NL + j] = 0;

			rtTmpAccess(F_OFFSET + i * NL + j);
			
            for (k = 0; k < NM; ++k)
            {
                F[i * NL + j] += C[i * NM + k] * D[k * NL + j];
				rtTmpAccess(C_OFFSET + i * NM + k);
				rtTmpAccess(D_OFFSET + k * NL + j);
				rtTmpAccess(F_OFFSET + i * NL + j);
				rtTmpAccess(F_OFFSET + i * NL + j);
			}
        }
    }

    /* G := E*F */
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NL; j++)
        {
            G[i * NL + j] = 0;

			rtTmpAccess(G_OFFSET + i * NL + j);

            for (k = 0; k < NJ; ++k)
            {
                G[i * NL + j] += E[i * NJ + k] * F[k * NL + j];
				rtTmpAccess(E_OFFSET + i * NJ + k);
				rtTmpAccess(F_OFFSET + k * NL + j);
				rtTmpAccess(G_OFFSET + i * NL + j);
				rtTmpAccess(G_OFFSET + i * NL + j);
			}
        }
    }

    return;
}



int main() {
    int* a = (int*)malloc( NK * NI * sizeof(int) );
    int* b = (int*)malloc( NK * NJ * sizeof(int) );
    int* c = (int*)malloc( NJ * NM * sizeof(int) );
    int* d = (int*)malloc( NM * NL * sizeof(int) );
    int* e = (int*)malloc( NI * NJ * sizeof(int) );
    int* f = (int*)malloc( NJ * NL * sizeof(int) );
    int* g = (int*)malloc( NI * NL * sizeof(int) );

    for (int i = 0; i < NK * NI; i++) {
            a[i] = i % 256;
    }

	for (int i = 0; i < NK * NJ; i++) {
            b[i] = i % 256;
    }

	for (int i = 0; i < NJ * NM; i++) {
            c[i] = i % 256;
    }

	for (int i = 0; i < NM * NL; i++) {
            d[i] = i % 256;
    }

#ifdef RD
    InitRD();
#endif
    
    mm3_cpu_trace(a, b, c, d, e, f, g);
	
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

