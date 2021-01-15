#include "data_size.h"

#ifdef PROFILE_RT
    #include "rt.h"
#endif

#ifdef RD
    #include "reda-spatial.h"
#endif

#if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET) && !defined(MEDIUM_DATASET)
    #define LARGE_DATASET
#endif
#ifdef MINI_DATASET
    #define NI 32
    #define NJ 32
    #define NK 32
    #define NL 32
    #define NM 32
#endif
#ifdef SMALL_DATASET
    #define NI 128
    #define NJ 128
    #define NL 128
    #define NK 128
    #define NM 128
#endif
#ifdef MEDIUM_DATASET
    #define NI 1024
    #define NJ 1024
    #define NK 1024
    #define NL 1024
    #define NM 1024
#endif
#ifdef LARGE_DATASET
    #define NI 2048
    #define NJ 2048
    #define NL 2048
    #define NK 2048
    #define NM 2048
#endif
#ifdef EXTRALARGE_DATASET
    #define NI 4096
    #define NJ 4096
    #define NL 4096
    #define NK 4096
    #define NM 4096
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
#ifdef PAPI_TIMER
    // Get starting timepoint
    PAPI_timer_init();
    PAPI_timer_start();
#endif
    
    mm3_cpu_trace(a, b, c, d, e, f, g);
	
#ifdef PROFILE_RT
    RTtoMR_AET();
#endif 
#ifdef PAPI_TIMER
    // Get ending timepoint
    PAPI_timer_end();
    PAPI_timer_print();
#endif
#ifdef PROFILE_RT
    dumpRtTmp();
    dumpMR();
#endif
    
#ifdef RD
    FiniRD();
#endif

    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);
    free(g);

    return 0;
}

