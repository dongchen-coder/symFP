
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"

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
            rtTmpAccess(V_OFFSET + 0 * N + i, 0, 0);
            rtTmpAccess(P_OFFSET + i * N + 0, 1, 1);
            rtTmpAccess(V_OFFSET + 0 * N + i, 2, 0);
            rtTmpAccess(Q_OFFSET + i * N + 0, 3, 2);
            for (j=1; j< N-1; j++) {
                p[i * N + j] = -c / (a*p[i * N + j-1]+b);
                q[i * N + j] = (-d*u[j * N + i-1]+(1.0 + 2.0 * d)*u[j * N + i] - f*u[j * N + i+1]-a*q[i * N + j-1])/(a*p[i * N + j-1]+b);
                rtTmpAccess(P_OFFSET + i * N + j-1, 4, 1);
                rtTmpAccess(P_OFFSET + i * N + j, 5, 1);
                rtTmpAccess(U_OFFSET + j * N + i-1, 6, 3);
                rtTmpAccess(U_OFFSET + j * N + i, 7, 3);
                rtTmpAccess(U_OFFSET + j * N + i+1, 8, 3);
                rtTmpAccess(Q_OFFSET + i * N + j-1, 9, 2);
                rtTmpAccess(P_OFFSET + i * N + j-1, 10, 1);
                rtTmpAccess(Q_OFFSET + i * N + j, 11, 2);
            }
            v[(N - 1) * N + i] = 1.0;
            rtTmpAccess(V_OFFSET + (N - 1) * N + i, 12, 0);
            for (j= N-2; j>=1; j--) {
                v[j * N + i] = p[i * N + j] * v[(j+1) * N + i] + q[i * N + j];
                rtTmpAccess(P_OFFSET + i * N + j, 13, 1);
                rtTmpAccess(V_OFFSET + (j+1) * N + i, 14, 0);
                rtTmpAccess(Q_OFFSET + i * N + j, 15, 2);
                rtTmpAccess(V_OFFSET + j * N + i, 16, 0);
            }
        }
        //Row Sweep
        for (i=1; i < N - 1; i++) {
            u[i * N + 0] = 1.0;
            p[i * N + 0] = 0.0;
            q[i * N + 0] = u[i * N + 0];
            rtTmpAccess(U_OFFSET + i * N + 0, 17, 3);
            rtTmpAccess(P_OFFSET + i * N + 0, 18, 1);
            rtTmpAccess(U_OFFSET + i * N + 0, 19, 3);
            rtTmpAccess(Q_OFFSET + i * N + 0, 20, 2);
            for (j=1; j< N - 1; j++) {
                p[i * N + j] = -f / (d*p[i * N + j-1]+e);
                q[i * N + j] = (-a*v[(i-1) * N + j]+(1.0 + 2.0 * a)*v[i * N + j] - c*v[(i+1) * N + j]-d*q[i * N + j-1])/(d*p[i * N + j-1]+e);
                rtTmpAccess(P_OFFSET + i * N + j-1, 21, 1);
                rtTmpAccess(P_OFFSET + i * N + j, 22, 1);
                rtTmpAccess(V_OFFSET + (i-1) * N + j, 23, 0);
                rtTmpAccess(V_OFFSET + i * N + j, 24, 0);
                rtTmpAccess(V_OFFSET + (i+1) * N + j, 25, 0);
                rtTmpAccess(Q_OFFSET + i * N + j-1, 26, 2);
                rtTmpAccess(P_OFFSET + i * N + j-1, 27, 1);
                rtTmpAccess(Q_OFFSET + i * N + j, 28, 2);
            }
            u[i * N + N - 1 ] = 1.0;
            rtTmpAccess(U_OFFSET + i * N + N - 1, 29, 3);
            for (j= N - 2; j>=1; j--) {
                u[i * N + j] = p[i * N + j] * u[i * N + j+1] + q[i * N + j];
                rtTmpAccess(P_OFFSET + i * N + j, 30, 1);
                rtTmpAccess(U_OFFSET + i * N + j+1, 31, 3);
                rtTmpAccess(Q_OFFSET + i * N + j, 32, 2);
                rtTmpAccess(U_OFFSET + i * N + j, 33, 3);
            }
        }
    }
    
}

int main() {

	double * p = (double *) malloc(N * N * sizeof(double));
	double * q = (double *) malloc(N * N * sizeof(double));
	double * v = (double *) malloc(N * N * sizeof(double));
	double * u = (double *) malloc(N * N * sizeof(double));

	adi_trace(p, q, v, u);
	
	OSL_ref(0);

	return 0;
}


