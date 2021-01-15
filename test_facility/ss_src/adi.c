#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif
#ifdef MINI_DATASET
	#define TSTEPS  2
	#define N       32
#endif
#ifdef SMALL_DATASET
	#define TSTEPS  10
    #define N       1024
#endif 
#ifdef MEDIUM_DATASET
	#define TSTEPS  50
	#define N       2048
#endif
#ifdef LARGE_DATASET
	#define TSTEPS  50
	#define N       8192
#endif
#ifdef EXTRALARGE_DATASET
	#define TSTEPS  100
	#define N       100000
#endif

#else
	#define TSTEPS  1
	#define N       4
#endif
void adi(double* p, double* q, double* v, double* u) {

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

	//Column Sweep
	for (i=0; i< N-2; i++) {
		v[0 * N + i+1] = 1.0;
		p[(i+1) * N + 0] = 0.0;
		q[(i+1) * N + 0] = v[0 * N + (i+1)];
		for (j=0; j< N-2; j++) {
			p[(i+1) * N + j+1] = -c / (a*p[(i+1) * N + j]+b);
			q[(i+1) * N + j+1] = (-d*u[(j+1) * N + i] + (1.0 + 2.0 * d) * u[(j+1) * N + i+1] - f * u[(j+1) * N + i+2] - a * q[(i+1) * N + j]) / (a*p[(i+1) * N + j]+b);
		}
		v[(N - 1) * N + (i+1)] = 1.0;
		for (j = 0; j < N-2; j++) {
			v[(N-2-j) * N + (i+1)] = p[(i+1) * N + (N-2-j)] * v[(N-1-j) * N + i+1] + q[(i+1) * N + N-2-j];
		}
	}
	//Row Sweep
	for (i=0; i < N - 2; i++) {
		u[(i+1) * N + 0] = 1.0;
		p[(i+1) * N + 0] = 0.0;
		q[(i+1) * N + 0] = u[(i+1) * N + 0];
		for (j=0; j< N - 2; j++) {
			p[(i+1) * N + j+1] = -f / (d*p[(i+1) * N + j]+e);
			q[(i+1) * N + j+1] = (-a * v[i * N + j+1] + (1.0 + 2.0 * a) * v[(i+1) * N + j+1] - c * v [(i+2) * N + j+1]- d * q[(i+1) * N + j]) / (d * p[(i+1) * N + j]+e);
		}
		u[(i+1) * N + N - 1 ] = 1.0;
		for (j=0; j<N-2; j++) {
			u[(i+1) * N + N-2-j] = p[(i+1) * N + N-2-j] * u[(i+1) * N + N-1-j] + q[(i+1) * N + N-2-j];
		}
	}

}
