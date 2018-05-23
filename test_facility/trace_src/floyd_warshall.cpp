#include "../utility/rt.h"
#include "../utility/data_size.h"

#ifdef ORG
	#define N 1024
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define PATH_OFFSET 0

void floyd_warshall_trace(double* path) {

	int k, i, j;

	for (k = 0; k < N; k++) 
		for(i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				//path[i * N + j] = path[i * N + j] < path[i * N + k] + path[k * N + j] ? path[i * N + j] : path[i * N + k] + path[k * N + j];
				rtTmpAccess(PATH_OFFSET + i * N + j);
				rtTmpAccess(PATH_OFFSET + i * N + k);
				rtTmpAccess(PATH_OFFSET + k * N + j);
				if (path[i * N + j] < path[i * N + k] + path[k * N + j]) {
					path[i * N + j] = path[i * N + j];
					rtTmpAccess(PATH_OFFSET + i * N + j);
					rtTmpAccess(PATH_OFFSET + i * N + j);
				} else {
					path[i * N + j] = path[i * N + k] + path[k * N + j];
					rtTmpAccess(PATH_OFFSET + i * N + k);
					rtTmpAccess(PATH_OFFSET + k * N + j);
					rtTmpAccess(PATH_OFFSET + i * N + j);
				}
			}
}

int main() {

	double* path = (double *)malloc(N * N * sizeof(double));

	floyd_warshall_trace(path);

	dumpRtTmp();
    RTtoMR_AET();
    dumpMR();

	return 0;
}

