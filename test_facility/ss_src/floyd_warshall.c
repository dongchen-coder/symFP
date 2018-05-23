
#define N 1024

void floyd_warshall(double* path) {

	int k, i, j;

	for (k = 0; k < N; k++) {
		for(i = 0; i < N; i++)
			for (j = 0; j < N; j++)
				path[i * N + j] = path[i * N + j] < path[i * N + k] + path[k * N + j] ? path[i * N + j] : path[i * N + k] + path[k * N + j];
	}

}
