
#define NR 256
#define NQ 256
#define NP 256

void doitgen(double* sum, double* A, double* C4) {

	int r, q, p, s;

	for (r = 0; r < NR; r++) {
		for (q = 0; q < NQ; q++)  {
			for (p = 0; p < NP; p++)  {
				sum[p] = 0.0;
				for (s = 0; s < NP; s++) {
					sum[p] += A[r * NQ * NP + q * NP + s] * C4[s * NP + p];
				}
			}
			for (p = 0; p < NP; p++) {
				A[r * NQ * NP + q * NP + p] = sum[p];
			}
    	}
	}

}

