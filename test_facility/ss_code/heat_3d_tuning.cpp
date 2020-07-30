
 /* Start to analysis array index
Array index info: Total number of references: 22
A.addr (((((i * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
A.addr ((((((i + 2) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
B.addr (((((i * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 2)
A.addr ((((((i + 1) * 2048) * 2048) + ((j + 2) * 2048)) + k) + 1)
B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
A.addr ((((((i + 1) * 2048) * 2048) + (j * 2048)) + k) + 1)
B.addr ((((((i + 1) * 2048) * 2048) + ((j + 2) * 2048)) + k) + 1)
B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
B.addr ((((((i + 1) * 2048) * 2048) + (j * 2048)) + k) + 1)
A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
A.addr (((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k)
A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
B.addr ((((((i + 2) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 2)
B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
B.addr (((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k)
B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %B
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 2046)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 2046)
----Loop inc: (j + 1)
----Loop predicate: <
------k
------Loop Bound: (0, 2046)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((((((i + 2) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access A.addr (((((i * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access A.addr ((((((i + 1) * 2048) * 2048) + ((j + 2) * 2048)) + k) + 1)
--------array access A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access A.addr ((((((i + 1) * 2048) * 2048) + (j * 2048)) + k) + 1)
--------array access A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 2)
--------array access A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access A.addr (((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k)
--------array access A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--i
--Loop Bound: (0, 2046)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 2046)
----Loop inc: (j + 1)
----Loop predicate: <
------k
------Loop Bound: (0, 2046)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access B.addr ((((((i + 2) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access B.addr (((((i * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access B.addr ((((((i + 1) * 2048) * 2048) + ((j + 2) * 2048)) + k) + 1)
--------array access B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access B.addr ((((((i + 1) * 2048) * 2048) + (j * 2048)) + k) + 1)
--------array access B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 2)
--------array access B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access B.addr (((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k)
--------array access B.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)
--------array access A.addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)

Finish analysis loops */ 
/* # of Out-most Loops: 2 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 20
------Sample number: 418
--------Sample number: 8564
----Sample number: 20
------Sample number: 418
--------Sample number: 8564
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* A_addr2	46047276 */
/* A_addr0	46047276 */
/* A_addr4	46047276 */
/* A_addr1	46047276 */
/* B_addr13	46047276 */
/* A_addr6	46047276 */
/* A_addr3	46047276 */
/* B_addr10	46047276 */
/* A_addr5	46047276 */
/* B_addr14	46047276 */
/* B_addr15	46047276 */
/* B_addr16	46047276 */
/* A_addr7	46047276 */
/* A_addr8	46047276 */
/* A_addr9	46047276 */
/* A_addr21	46047276 */
/* B_addr11	46047276 */
/* B_addr12	46047276 */
/* B_addr17	46047276 */
/* B_addr18	46047276 */
/* B_addr19	46047276 */
/* B_addr20	46047276 */
#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#ifdef PAPI_TIMER
#  include "papi_timer.h"
#endif
using namespace std;
using namespace placeholders;
map<uint64_t, double> RT;
map<uint64_t, double> MR;
uint64_t sample_sum = 0;
double total_reuse, total_neighbor, total_scale, total_fold, total_interchunk;
void rtHistoCal( map<uint64_t, double> &rth, uint64_t rt, double val ) {
    if (rth.find(rt) == rth.end()) { 
        rth[rt] = val;
    } else {
        rth[rt] += val;
    }
    return;
}
void subBlkRT(map<uint64_t, double> &rth, int rt, double cnt) {
    int msb = 0;
    int tmp_rt = rt;
    while(tmp_rt != 0) {
        tmp_rt = tmp_rt / 2;
        ++msb;
    }
    if (msb >= BIN_SIZE) {
        int diff = (pow(2, msb) - pow(2, msb-1)) / BIN_SIZE;
        for (int b = pow(2, msb-1); b <= pow(2, msb); b+=diff) {
            if (rt < b) {
                rtHistoCal(rth, b - diff, cnt);
                break;
            }
        }
    }
    else {
        rtHistoCal(rth, pow(2, msb-1), cnt);
    }
    return;
}
int getChunkNum(uint64_t is, int thread_cnt) {
    return is / (CHUNK_SIZE * thread_cnt);
}
int getChunkID(uint64_t i, int thread_cnt) {
    return floor(i / (CHUNK_SIZE * thread_cnt));
}
int getThreadID(uint64_t i, int thread_cnt) {
    return i / CHUNK_SIZE - floor(i / (CHUNK_SIZE * thread_cnt))*thread_cnt ;
}
int getThreadLocalPos(uint64_t i) {
    return i % CHUNK_SIZE;
}
int search_src_candidate_neighbor(uint64_t i, int thread_cnt, function<uint64_t(uint64_t)> calAddr) {
    int c_start = i % CHUNK_SIZE + i / (CHUNK_SIZE * thread_cnt) * thread_cnt * CHUNK_SIZE;
    int c_end = c_start + (thread_cnt - 1) * CHUNK_SIZE;
    for (int c = i + CHUNK_SIZE; c <= c_end; c=c+CHUNK_SIZE) {
        if (calAddr(i) == calAddr(c)) { return getThreadID(c, thread_cnt); }
    }
    return -1;
}
int search_sink_candidate_neighbor(uint64_t i, int thread_cnt, function<uint64_t(uint64_t)> calAddr) {
    int c_start = i % CHUNK_SIZE + i / (CHUNK_SIZE * thread_cnt) * thread_cnt * CHUNK_SIZE;
    for (int c = c_start; c <= i; c=c+CHUNK_SIZE) {
        if (calAddr(i) == calAddr(c)) { return getThreadID(c, thread_cnt); }
    }
    return -1;
}
uint64_t parallel_predict(uint64_t i_src, uint64_t i_sink, uint64_t rt, uint64_t lsrc, uint64_t lsink, uint64_t middle_accesses, int thread_cnt, bool is_normal_ref, bool is_in_same_loop, int & type, function<uint64_t(uint64_t)> srcAddrCal, function<uint64_t(uint64_t)> sinkAddrCal) {
    total_reuse += 1.0;
    uint64_t parallel_rt = rt;
    int tsrc = getThreadID(i_src, thread_cnt);
    int tsink = getThreadID(i_sink, thread_cnt);
    int dT = tsink - tsrc;
    if (!is_in_same_loop || getChunkID(i_src, thread_cnt) != getChunkID(i_sink, thread_cnt)) {
#ifdef DEBUG
        cout << "Inter Chunk Reuse" << endl;
        ;
        cout << "rt " << rt << endl;
#endif
    type = 0; // code for inter chunk reuse
        parallel_rt = rt * thread_cnt - CHUNK_SIZE * thread_cnt * (lsrc*(thread_cnt - tsrc) + lsink * tsink) + CHUNK_SIZE * thread_cnt * lsrc - (thread_cnt - 1) * middle_accesses + dT;
    } else if (!is_normal_ref) {
        /* intra chunk reuse */
#ifdef DEBUG
        cout << "Neighboring Effect" << endl;
#endif
        total_neighbor += 1.0;
        int tsrc_neighbor = search_src_candidate_neighbor(i_src, thread_cnt, srcAddrCal);
        int tsink_neighbor = search_sink_candidate_neighbor(i_sink, thread_cnt, sinkAddrCal);
        if (tsrc_neighbor >= 0) {
#ifdef DEBUG
            cout << "Find sink in src neighbor at " << tsrc_neighbor << endl;
            cout << "Neighbor Effect: " << tsrc_neighbor - tsrc << endl;
#endif
        type = 1; // code for src neighboring effect
            return tsrc_neighbor - tsrc;
        } else if (tsink_neighbor >= 0) {
#ifdef DEBUG
            cout << "Find sink in sink neighbor at " << tsink_neighbor << endl;
            cout << "Neighbor Effect: " << rt * thread_cnt + tsink_neighbor - tsink << endl;
#endif
            if (getChunkID(i_src, thread_cnt) == getChunkID(i_sink, thread_cnt)) {
        type = 2; // code for sink neighboring effect
                return rt * thread_cnt + tsink_neighbor - tsink;
            }
        }
    } else if (getChunkID(i_src, thread_cnt) == getChunkID(i_sink, thread_cnt)) {
        /* same thread -- scaling effect */
        if (dT == 0) {
#ifdef DEBUG
            cout << "Scaling Effect" << endl;
#endif
            total_scale += 1.0;
            parallel_rt = rt * thread_cnt;
            type = 3; // code for scaling effect
        } else if (getThreadLocalPos(i_src) <= getThreadLocalPos(i_sink)) { // src-sink order
            if ((rt * thread_cnt - CHUNK_SIZE * lsrc * thread_cnt * dT + dT) < 0) { printf("NORMAL ORDER NEGATIVE PRI\n"); }
#ifdef DEBUG
            cout << "Src-Sink Order Folding Effect" << endl;
#endif
            type = 4; // code for src-sink order folding effect
            total_fold += 1.0;
            parallel_rt = rt * thread_cnt - CHUNK_SIZE * lsrc * thread_cnt * dT + abs(dT);
        } else { // sink-src order
            if ((rt * thread_cnt - CHUNK_SIZE * lsrc * thread_cnt * dT + dT) < 0) { printf("REVERSE ORDER NEGATIVE PRI\n"); }
#ifdef DEBUG
            cout << "Sink-Src Order Folding Effect" << endl;
#endif
            // parallel_rt = CHUNK_SIZE * lsrc * thread_cnt * dT - (rt * thread_cnt) - abs(dT);
            type = 5; // code for sink-src order folding effect
            return 0;
        }
    }
    return parallel_rt;
}
void RTtoMR_AET() {
     map<uint64_t, double> P;
    double total_num_RT = 0;
    uint64_t max_RT = 0;
    for ( map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
        total_num_RT += it->second;
        if (max_RT < it->first) {
            max_RT = it->first;
        }
    }
    double accumulate_num_RT = 0;
    for ( map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
        P[it->first] = accumulate_num_RT / total_num_RT;
        accumulate_num_RT += it->second;
    }
    P[0] = 1;
    double sum_P = 0;
    uint64_t t = 0;
    uint64_t prev_t = 0;
    for (uint64_t c = 0; c <= max_RT && c <= 327680; c++) {
        while (sum_P < c && t <= max_RT) {
            if (P.find(t) != P.end()) {
                sum_P += P[t];
                prev_t = t;
            } else {
                sum_P += P[prev_t];
            }
            t++;
        }
        MR[c] = P[prev_t];
    }
    return;
}
void rtDump() {
    cout << "Start to dump reuse time histogram\n";
    for (map<uint64_t, double>::iterator it = RT.begin(), eit = RT.end(); it != eit; ++it) {
        cout << it->first << ", " << it->second << "\n";
    }
    return;
}
void dumpMR() {
    cout << "miss ratio" << endl;
     map<uint64_t, double>::iterator it1 = MR.begin();
     map<uint64_t, double>::iterator it2 = MR.begin();
    while(it1 != MR.end()) {
        while(1) {
             map<uint64_t, double>::iterator it3 = it2;
            ++it3;
            if (it3 == MR.end()) {
                break;
            }
            if (it1->second - it3->second < 0.00001) {
                ++it2;
            } else {
                break;
            }
        }
        cout << it1->first << ", " << it1->second << endl;
        if (it1 != it2) {
            cout << it2->first << ", " << it2->second << endl;
        }
        it1 = ++it2;
        it2 = it1;
    }
    return;
}
/* Dump the reuse statistics */
void statDump() {
    cout << "Total Neighboring Reuses: " << total_neighbor / total_reuse << endl;
    cout << "Total Scaling Reuses: " << total_scale / total_reuse << endl;
    cout << "Total Folding Src-Sink Reuses: " << total_fold / total_reuse << endl;
    cout << "Total Inter Chunk Reuses: " << total_interchunk / total_reuse << endl;
    cout << "Total Reuses: " << total_reuse << endl;
    return;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr ((((((i + 2) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 0 */
int calAddrA_addr0( int i, int j, int k) {
    int result = (((((((i + 2) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 1 */
int calAddrA_addr1( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array A_addr	j i k */ 
/* i */
/* i */
/* i, j */
/* A_addr (((((i * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 2 */
int calAddrA_addr2( int i, int j, int k) {
    int result = ((((((i * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr ((((((i + 1) * 2048) * 2048) + ((j + 2) * 2048)) + k) + 1) 3 */
int calAddrA_addr3( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 2) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 4 */
int calAddrA_addr4( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr ((((((i + 1) * 2048) * 2048) + (j * 2048)) + k) + 1) 5 */
int calAddrA_addr5( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + (j * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 2) 6 */
int calAddrA_addr6( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 2)) * 8 / 64;
    return result;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 7 */
int calAddrA_addr7( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr (((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) 8 */
int calAddrA_addr8( int i, int j, int k) {
    int result = ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k)) * 8 / 64;
    return result;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 9 */
int calAddrA_addr9( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 10 */
int calAddrB_addr10( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr ((((((i + 2) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 11 */
int calAddrB_addr11( int i, int j, int k) {
    int result = (((((((i + 2) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 12 */
int calAddrB_addr12( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array B_addr	j i k */ 
/* i */
/* i */
/* i, j */
/* B_addr (((((i * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 13 */
int calAddrB_addr13( int i, int j, int k) {
    int result = ((((((i * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr ((((((i + 1) * 2048) * 2048) + ((j + 2) * 2048)) + k) + 1) 14 */
int calAddrB_addr14( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 2) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 15 */
int calAddrB_addr15( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr ((((((i + 1) * 2048) * 2048) + (j * 2048)) + k) + 1) 16 */
int calAddrB_addr16( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + (j * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 2) 17 */
int calAddrB_addr17( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 2)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 18 */
int calAddrB_addr18( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr (((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) 19 */
int calAddrB_addr19( int i, int j, int k) {
    int result = ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k)) * 8 / 64;
    return result;
}
/* Array B_addr	i j k */ 
/* i */
/* B_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 20 */
int calAddrB_addr20( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
/* Array A_addr	i j k */ 
/* i */
/* A_addr ((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1) 21 */
int calAddrA_addr21( int i, int j, int k) {
    int result = (((((((i + 1) * 2048) * 2048) + ((j + 1) * 2048)) + k) + 1)) * 8 / 64;
    return result;
}
void ref_A_addr2(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr2 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr0(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr0 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr4(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr4 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr1(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr1;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr13(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr13 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr6(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr6 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr3(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr3 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr10(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr10( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr10, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr10;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr10(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr10 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr10;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr5(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr5(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr5 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr14(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr14 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr15(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr15(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr15, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr15 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr15;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr16(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr16(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr16, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr16 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr16;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr7(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr7 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr8(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr8 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr9(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 2046; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 2046; j=(j + 1)) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr1;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr2, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr4, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr6, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr9, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = (getChunkNum(2046, thread_cnt) - getChunkID((i_Start -0), thread_cnt)) * CHUNK_SIZE * thread_cnt * 46047276 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 46047276 + 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr9 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr21(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr21( i, j, k) == calAddrA_addr21(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr21, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr21, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[A_addr21 --> A_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "A_addr21;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "A_addr21;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr11(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr11 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr12(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr12 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr17(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr17(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr17, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr17 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr17;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr18(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr18(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr18, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr18 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr18;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr19(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr19(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr19, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr19 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr19;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr20(int thread_cnt) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 8564;) {
SAMPLE:
        int i_Start = rand() % (2046 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 2046) { goto SAMPLE; }
        if ( (2046 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (2046 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (2046 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (2046 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 2046; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 2046; j=(j + 1)) {
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 2046; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr11( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr11, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr12( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr12, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr13( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr13, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr14( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr14, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr15( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr15, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr15;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr16( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr16, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr16;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr17( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr17, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr17;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr18( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr18, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr18;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr19( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr19, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr19;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr20( i, j, k) == calAddrB_addr20(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr20, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr20, _1, j, k);
                                /* i 46047276 */
                                /* j 22506 */
                                /* k 11 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 46047276, 46047276, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr20 --> B_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr20;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr20;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
int main() {
#ifdef PAPI_TIMER
    PAPI_timer_init();
    PAPI_timer_start();
#endif
    int tlb = THREAD_LB;
    int tub = THREAD_UB;
    /* metadata to derive best thread */
    double accumulate_num_RT, expected_rt;
    double min_expected_rt = 0.0;
    int min_tnum = tlb;
    for (int t = tlb; t <= tub; t++) {
        accumulate_num_RT = 0.0;
        accumulate_num_RT = 0.0;
        total_reuse = 0.0;
        total_neighbor = 0.0;
        total_scale = 0.0;
        total_fold = 0.0;
        total_interchunk = 0.0;
        expected_rt = 0.0;
        /* Currently we consider all even number threads only */
        if ( t % 2 != 0 ) { continue; }        ref_A_addr2(t);
        ref_A_addr0(t);
        ref_A_addr4(t);
        ref_A_addr1(t);
        ref_B_addr13(t);
        ref_A_addr6(t);
        ref_A_addr3(t);
        ref_B_addr10(t);
        ref_A_addr5(t);
        ref_B_addr14(t);
        ref_B_addr15(t);
        ref_B_addr16(t);
        ref_A_addr7(t);
        ref_A_addr8(t);
        ref_A_addr9(t);
        ref_A_addr21(t);
        ref_B_addr11(t);
        ref_B_addr12(t);
        ref_B_addr17(t);
        ref_B_addr18(t);
        ref_B_addr19(t);
        ref_B_addr20(t);
        /* Compute expected reuse for current thread number */
        for ( map<uint64_t, double>::iterator it = RT.begin(), eit = RT.end(); it != eit; ++it) {
            accumulate_num_RT += it->second;
        }
        for ( map<uint64_t, double>::iterator it = RT.begin(), eit = RT.end(); it != eit; ++it) {
            expected_rt += (it->first * it->second / accumulate_num_RT);
        }
        /* assign min_tnum if its expected reuse is shorter */
        if (t == tlb)
            min_expected_rt = expected_rt;
        else if (expected_rt <= min_expected_rt)
            min_tnum = t;
        /* remove all RT and collect it for another thread id */
        statDump();
        RT.clear();
    }
    cout << "best thread number is: " << min_tnum << endl;
#ifdef PAPI_TIMER
    PAPI_timer_end();
    PAPI_timer_print();
#endif
    return 0;
}
 /* Analyze function: heat_3d */ 
