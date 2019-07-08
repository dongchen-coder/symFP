
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
				rtTmpAccess(PATH_OFFSET + i * N + j, 0, 0);
				rtTmpAccess(PATH_OFFSET + i * N + k, 1, 0);
				rtTmpAccess(PATH_OFFSET + k * N + j, 2, 0);
				if (path[i * N + j] < path[i * N + k] + path[k * N + j]) {
					path[i * N + j] = path[i * N + j];
					rtTmpAccess(PATH_OFFSET + i * N + j, 3, 0);
					rtTmpAccess(PATH_OFFSET + i * N + j, 4, 0);
				} else {
					path[i * N + j] = path[i * N + k] + path[k * N + j];
					rtTmpAccess(PATH_OFFSET + i * N + k, 5, 0);
					rtTmpAccess(PATH_OFFSET + k * N + j, 6, 0);
					rtTmpAccess(PATH_OFFSET + i * N + j, 7, 0);
				}
			}
}

int main() {

	double* path = (double *)malloc(N * N * sizeof(double));

	floyd_warshall_trace(path);

	dumpSetSize();

	return 0;
}
