#include "../utility/data_size.h"

#ifdef PROFILE_RT
#include "../utility/rt.h"
#endif

#ifdef RD
#include "../utility/reda-spatial.h"
#endif

#ifdef ORG
	#define N 1024
	#define TSTEPS 10
#elif defined (TX)
#elif defined (FX)
#elif defined (EX)
#endif

#define P_OFFSET 0
#define Q_OFFSET 0 + N * N
#define V_OFFSET 0 + N * N + N * N 
#define U_OFFSET 0 + N * N + N * N + N * N

void adi_trace(double* p, double* q, double* v, double* u) {

	int t, i, j;
	double DX, DY, DT, B1, B2, mul1, mul2, a, b, c, d, e, f;

	DX = 1.0 / N;
	DY = 1.0 / N;
	DT = 1.0 / TSTEPS;
	B1 = 2.0;
	B2 = 1.0;
	mul1 = B1 * DT / (DX * DX);
	mul2 = B2 * DT / (DY * DY);

	a = -mul1 / 2.0;
	b = 1.0 + mul1;
	c = a;
	d = -mul2 / 2.0;
	e = 1.0 + mul2;
	f = d;

	for (t=1; t <= TSTEPS; t++) {
	//Column Sweep
		for (i=1; i< N-1; i++) {
			v[0 * N + i] = 1.0;
			p[i * N + 0] = 0.0;
			q[i * N + 0] = v[0 * N + i];
			rtTmpAccess(V_OFFSET + 0 * N + i);
			rtTmpAccess(P_OFFSET + i * N + 0);
			rtTmpAccess(V_OFFSET + 0 * N + i);
			rtTmpAccess(Q_OFFSET + i * N + 0);
			for (j=1; j< N-1; j++) {
				p[i * N + j] = -c / (a*p[i * N + j-1]+b);
				q[i * N + j] = (-d*u[j * N + i-1]+(1.0 + 2.0 * d)*u[j * N + i] - f*u[j * N + i+1]-a*q[i * N + j-1])/(a*p[i * N + j-1]+b);
				rtTmpAccess(P_OFFSET + i * N + j-1);
				rtTmpAccess(P_OFFSET + i * N + j);
				rtTmpAccess(U_OFFSET + j * N + i-1);
				rtTmpAccess(U_OFFSET + j * N + i);
				rtTmpAccess(U_OFFSET + j * N + i+1);
				rtTmpAccess(Q_OFFSET + i * N + j-1);
				rtTmpAccess(P_OFFSET + i * N + j-1);
				rtTmpAccess(Q_OFFSET + i * N + j);
			}
			v[(N - 1) * N + i] = 1.0;
			rtTmpAccess(V_OFFSET + (N - 1) * N + i);
			for (j= N-2; j>=1; j--) {
				v[j * N + i] = p[i * N + j] * v[(j+1) * N + i] + q[i * N + j];
				rtTmpAccess(P_OFFSET + i * N + j);
				rtTmpAccess(V_OFFSET + (j+1) * N + i);
				rtTmpAccess(Q_OFFSET + i * N + j);
				rtTmpAccess(V_OFFSET + j * N + i);
			}
		}
		//Row Sweep
		for (i=1; i < N - 1; i++) {
			u[i * N + 0] = 1.0;
			p[i * N + 0] = 0.0;
			q[i * N + 0] = u[i * N + 0];
			rtTmpAccess(U_OFFSET + i * N + 0);
			rtTmpAccess(P_OFFSET + i * N + 0);
			rtTmpAccess(U_OFFSET + i * N + 0);
			rtTmpAccess(Q_OFFSET + i * N + 0);
			for (j=1; j< N - 1; j++) {
				p[i * N + j] = -f / (d*p[i * N + j-1]+e);
				q[i * N + j] = (-a*v[(i-1) * N + j]+(1.0 + 2.0 * a)*v[i * N + j] - c*v[(i+1) * N + j]-d*q[i * N + j-1])/(d*p[i * N + j-1]+e);
				rtTmpAccess(P_OFFSET + i * N + j-1);
				rtTmpAccess(P_OFFSET + i * N + j);
				rtTmpAccess(V_OFFSET + (i-1) * N + j);
				rtTmpAccess(V_OFFSET + i * N + j);
				rtTmpAccess(V_OFFSET + (i+1) * N + j);
				rtTmpAccess(Q_OFFSET + i * N + j-1);
				rtTmpAccess(P_OFFSET + i * N + j-1);
				rtTmpAccess(Q_OFFSET + i * N + j);
			}
			u[i * N + N - 1 ] = 1.0;
			rtTmpAccess(U_OFFSET + i * N + N - 1);
			for (j= N - 2; j>=1; j--) {
				u[i * N + j] = p[i * N + j] * u[i * N + j+1] + q[i * N + j];
				rtTmpAccess(P_OFFSET + i * N + j);
				rtTmpAccess(U_OFFSET + i * N + j+1);
				rtTmpAccess(Q_OFFSET + i * N + j);
				rtTmpAccess(U_OFFSET + i * N + j);
			}
		}
	}

}

int main() {

	double * p = (double *) malloc(N * N * sizeof(double));
	double * q = (double *) malloc(N * N * sizeof(double));
	double * v = (double *) malloc(N * N * sizeof(double));
	double * u = (double *) malloc(N * N * sizeof(double));

#ifdef RD
    InitRD();
#endif
    
	adi_trace(p, q, v, u);
	
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



