
#define _PB_TMAX 10
#define _PB_NY 1024
#define _PB_NX 1024

void fdtd_2d(double* ey, double* ex, double* hz) {

	int t, i, j;
    double _fict_;

    for (j = 0; j < _PB_NY; j++)
        ey[0 * _PB_NY + j] = _fict_;

	for (i = 1; i < _PB_NX; i++)
		for (j = 0; j < _PB_NY; j++)
			ey[i * _PB_NY + j] = ey[i * _PB_NY + j] - 0.5 * (hz[i * _PB_NY + j] - hz[(i-1) * _PB_NY + j]);

	for (i = 0; i < _PB_NX; i++)
		for (j = 1; j < _PB_NY; j++)
			ex[i * _PB_NY + j] = ex[i * _PB_NY + j] - 0.5 * (hz[i * _PB_NY + j] - hz[i * _PB_NY + j - 1]);

	for (i = 0; i < _PB_NX - 1; i++)
		for (j = 0; j < _PB_NY - 1; j++)
			hz[i * _PB_NY + j] = hz[i * _PB_NY + j] - 0.7 * (ex[i * _PB_NY + j + 1] - ex[i * _PB_NY + j] + ey[(i+1) * _PB_NY + j] - ey[i * _PB_NY + j]);

}
