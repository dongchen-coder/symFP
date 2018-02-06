
#define N 1024


void nussinov(double* table, double* seq) {

	int i, j, k;

	for (i = N - 1; i >= 0; i--) {
		for (j=i+1; j < N; j++) {
			if (j-1>=0)
				table[i * N + j] = max_score(table[i * N + j], table[i * N + j-1]);
			if (i+1 < N)
				table[i * N +j] = max_score(table[i * N + j], table[(i+1) * N + j]);
			if (j-1>=0 && i+1 < N) {
				/* don't allow adjacent elements to bond */
				if (i<j-1)
					table[i * N + j] = max_score(table[i * N + j], table[(i+1) * N + j-1]+match(seq[i], seq[j]));
				else
					table[i * N + j] = max_score(table[i * N + j], table[(i+1) * N + j-1]);
			}
			for (k=i+1; k<j; k++) {
				table[i * N + j] = max_score(table[i * N + j], table[i * N + k] + table[(k+1) * N + j]);
			}
		}
	}

}



