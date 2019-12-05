
#include "../utility/reference_lease.h"
#include "../utility/data_size.h"
#ifdef PAPI_TIMER
#   include "../utility/papi_timer.h"
#endif

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
            if (j-1>=0) {
                table[i * N + j] = max_score(table[i * N + j], table[i * N + j-1]);
                rtTmpAccess(TABLE_OFFSET + i * N + j, 0, 0);
                rtTmpAccess(TABLE_OFFSET + i * N + j-1, 1, 0);
                //if (table[i * N + j] > table[i * N + j-1])
                //    rtTmpAccess(TABLE_OFFSET + i * N + j, 2, 0);
                //else
                //    rtTmpAccess(TABLE_OFFSET + i * N + j-1, 3, 0);
                rtTmpAccess(TABLE_OFFSET + i * N + j, 2, 0);
            }
            if (i+1 < N) {
                table[i * N +j] = max_score(table[i * N + j], table[(i+1) * N + j]);
                rtTmpAccess(TABLE_OFFSET + i * N + j, 3, 0);
                rtTmpAccess(TABLE_OFFSET + (i+1) * N + j, 4, 0);
                //if (table[i * N + j] >= table[(i+1) * N + j])
                //    rtTmpAccess(TABLE_OFFSET + i * N + j, 7, 0);
                //else
                //    rtTmpAccess(TABLE_OFFSET + i * N + j, 8, 0);
                rtTmpAccess(TABLE_OFFSET + i * N + j, 5, 0);
            }
            if (j-1>=0 && i+1 < N) {
                /* don't allow adjacent elements to bond */
                if (i<j-1) {
                    table[i * N + j] = max_score(table[i * N + j], table[(i+1) * N + j-1]+match(seq[i], seq[j]));
                    rtTmpAccess(TABLE_OFFSET + i * N + j, 6, 0);
                    rtTmpAccess(TABLE_OFFSET + (i+1) * N + j-1, 7, 0);
                    rtTmpAccess(SEQ_OFFSET + i, 8, 1);
                    rtTmpAccess(SEQ_OFFSET + j, 9, 1);
                    
                    //if (table[i * N + j] >= table[(i+1) * N + j-1]+match(seq[i], seq[j])) {
                    //    rtTmpAccess(TABLE_OFFSET + i * N + j, 14, 0);
                    //} else {
                    //    rtTmpAccess(TABLE_OFFSET + (i+1) * N + j-1, 15, 0);
                    //    rtTmpAccess(SEQ_OFFSET + i, 16, 1);
                    //    rtTmpAccess(SEQ_OFFSET + j, 17, 1);
                    //}
                    rtTmpAccess(TABLE_OFFSET + i * N + j, 10, 0);
                } else {
                    table[i * N + j] = max_score(table[i * N + j], table[(i+1) * N + j-1]);
                    rtTmpAccess(TABLE_OFFSET + i * N + j, 11, 0);
                    rtTmpAccess(TABLE_OFFSET + (i+1) * N + j-1, 12, 0);
                    //if (table[i * N + j] >= table[(i+1) * N + j-1]) {
                    //    rtTmpAccess(TABLE_OFFSET + i * N + j, 21, 0);
                    //} else {
                    //    rtTmpAccess(TABLE_OFFSET + i * N + j, 22, 0);
                    //}
                    rtTmpAccess(TABLE_OFFSET + i * N + j, 13, 0);
                }
            }
            for (k=i+1; k<j; k++) {
                table[i * N + j] = max_score(table[i * N + j], table[i * N + k] + table[(k+1) * N + j]);
                rtTmpAccess(TABLE_OFFSET + i * N + j, 14, 0);
                rtTmpAccess(TABLE_OFFSET + i * N + k, 15, 0);
                rtTmpAccess(TABLE_OFFSET + (k+1) * N + j, 16, 0);
                //if (table[i * N + j] >= table[i * N + k] + table[(k+1) * N + j]) {
                //    rtTmpAccess(TABLE_OFFSET + i * N + j, 27, 0);
                //} else {
                //    rtTmpAccess(TABLE_OFFSET + i * N + k, 28, 0);
                //    rtTmpAccess(TABLE_OFFSET + (k+1) * N + j, 29, 0);
                //}
                rtTmpAccess(TABLE_OFFSET + i * N + j, 17, 0);
            }
        }
    }
}

int main() {

	double* table = (double *)malloc(N * N * sizeof(double));
	double* seq = (double *)malloc(N * sizeof(double));

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			table[i * N + j] = (i * N + j) % 128;
		}
		seq[i] = i % 64;
	}

#ifdef PAPI_TIMER
    PAPI_timer_init();
    PAPI_timer_start();
#endif
	nussinov_trace(table, seq);

#ifdef PAPI_TIMER
    PAPI_timer_end();
    PAPI_timer_print();
#endif

    printf("CARL Lease Assignment:\n");
#ifdef PAPI_TIMER
    PAPI_timer_start();
#endif
    OSL_ref(0);
#ifdef PAPI_TIMER
    PAPI_timer_end();
    PAPI_timer_print();
#endif
	return 0;
}
