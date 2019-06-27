#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

#ifdef ORG
	#define NX 1024
	#define NY 1024
#elif defined(TX)
	#define NX 1024
	#define NY 2048
#elif defined(FX)
	#define NX 1024
	#define NY 4096
#elif defined(EX)
	#define NX 1024
	#define NY 8192
#endif


#define A_OFFSET 0
#define S_OFFSET NX * NY
#define Q_OFFSET NX * NY + NY
#define R_OFFSET NX * NY + NY + NX
#define P_OFFSET NX * NY + NY + NX + NX


void bicg_cpu_trace(double* A, double* r, double* s, double* p, double* q, unsigned int nx, int ny) {
 
	int i,j;

    for (i = 0; i < NY; i++)
    {
        s[i] = 0.0;
    	rtTmpAccess(S_OFFSET + i);
	}

    for (i = 0; i < NX; i++)
    {
        q[i] = 0.0;
		rtTmpAccess(Q_OFFSET + i);
        for (j = 0; j < NY; j++)
        {
            s[j] = s[j] + r[i] * A[i * NY + j];
            q[i] = q[i] + A[i * NY + j] * p[j];
        	rtTmpAccess(S_OFFSET + j);
			rtTmpAccess(R_OFFSET + i);
			rtTmpAccess(A_OFFSET + i * NY + j);
			rtTmpAccess(S_OFFSET + j);
			rtTmpAccess(Q_OFFSET + i);
			rtTmpAccess(A_OFFSET + i * NY + j);
			rtTmpAccess(P_OFFSET + j);
			rtTmpAccess(Q_OFFSET + i);
		}
    }

    return;
}
 


int main() {
    double* A = (double*)malloc( (NX*NY)*sizeof(double));
    double* r = (double*)malloc( NX*sizeof(double));
    double* s = (double*)malloc( NY*sizeof(double));
    double* q = (double*)malloc( NX*sizeof(double));
    double* p = (double*)malloc( NY*sizeof(double));

    for (int i = 0; i < NX; i++) {
        r[i] = i % 256;
    }

    for (int i = 0; i< NY; i++) {
        p[i] = i % 256;
    }

    for (int i = 0; i< NY*NX; i++) {
        A[i] = i % 128;
    }
    
#ifdef RD
    InitRD();
#endif

    bicg_cpu_trace(A, r, s, p ,q, NX, NY);
    
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

