#include "../utility/rt.h"
#include "../utility/data_size.h"

#include <math.h>

#ifdef ORG
	#define N 1024
#elif defined(TX)
#elif defined(FX)
#elif defined(EX)
#endif

#define match(b1, b2) (((b1)+(b2)) == 3 ? 1 : 0)
#define max_score(s1, s2) ((s1 >= s2) ? s1 : s2)


#define TABLE_OFFSET 0
#define SEQ_OFFSET N * N

void nussinov_trace(double* table, double* seq) {

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

int main() {

	double* table = (double *)malloc(N * N * sizeof(double));
	double* seq = (double *)malloc(N * sizeof(double));

	nussinov_trace(table, seq);

	dumpRtTmp();
    RTtoMR_AET();
    dumpMR();

	return 0;
}

