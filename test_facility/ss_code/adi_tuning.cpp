
 /* Start to analysis array index
Array index info: Total number of references: 34
v.addr ((0 + i) + 1)
p.addr (((i + 1) * 8192) + j)
u.addr (((j + 1) * 8192) + i)
u.addr ((((j + 1) * 8192) + i) + 1)
p.addr (((i + 1) * 8192) + 0)
v.addr (0 + (i + 1))
q.addr (((i + 1) * 8192) + 0)
q.addr ((((i + 1) * 8192) + j) + 1)
v.addr (67100672 + (i + 1))
p.addr (((i + 1) * 8192) + (8190 - j))
p.addr ((((i + 1) * 8192) + j) + 1)
v.addr (((8190 - j) * 8192) + (i + 1))
u.addr (((i + 1) * 8192) + 0)
p.addr (((i + 1) * 8192) + 0)
u.addr (((i + 1) * 8192) + 0)
q.addr (((i + 1) * 8192) + 0)
p.addr (((i + 1) * 8192) + j)
u.addr ((((j + 1) * 8192) + i) + 2)
q.addr (((i + 1) * 8192) + j)
p.addr (((i + 1) * 8192) + j)
q.addr ((((i + 1) * 8192) + j) + 1)
u.addr ((((i + 1) * 8192) + 8192) - 1)
p.addr (((((i + 1) * 8192) + 8192) - 2) - j)
u.addr (((((i + 1) * 8192) + 8192) - 1) - j)
q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
u.addr (((((i + 1) * 8192) + 8192) - 2) - j)
v.addr ((((8191 - j) * 8192) + i) + 1)
q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
p.addr ((((i + 1) * 8192) + j) + 1)
v.addr (((i * 8192) + j) + 1)
v.addr ((((i + 1) * 8192) + j) + 1)
v.addr ((((i + 2) * 8192) + j) + 1)
q.addr (((i + 1) * 8192) + j)
p.addr (((i + 1) * 8192) + j)
BC Array cost info: Total number of arrays: 4
p.addr 24
q.addr 20
u.addr 19
v.addr 19
BC Average cost: 2.050000e+01

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %p
double* %q
double* %v
double* %u

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 8190)
--Loop inc: (i + 1)
--Loop predicate: <
----array access v.addr ((0 + i) + 1)
----array access p.addr (((i + 1) * 8192) + 0)
----array access v.addr (0 + (i + 1))
----array access q.addr (((i + 1) * 8192) + 0)
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p.addr (((i + 1) * 8192) + j)
------array access p.addr ((((i + 1) * 8192) + j) + 1)
------array access u.addr (((j + 1) * 8192) + i)
------array access u.addr ((((j + 1) * 8192) + i) + 1)
------array access u.addr ((((j + 1) * 8192) + i) + 2)
------array access q.addr (((i + 1) * 8192) + j)
------array access p.addr (((i + 1) * 8192) + j)
------array access q.addr ((((i + 1) * 8192) + j) + 1)
----array access v.addr (67100672 + (i + 1))
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p.addr (((i + 1) * 8192) + (8190 - j))
------array access v.addr ((((8191 - j) * 8192) + i) + 1)
------array access q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
------array access v.addr (((8190 - j) * 8192) + (i + 1))
--i
--Loop Bound: (0, 8190)
--Loop inc: (i + 1)
--Loop predicate: <
----array access u.addr (((i + 1) * 8192) + 0)
----array access p.addr (((i + 1) * 8192) + 0)
----array access u.addr (((i + 1) * 8192) + 0)
----array access q.addr (((i + 1) * 8192) + 0)
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p.addr (((i + 1) * 8192) + j)
------array access p.addr ((((i + 1) * 8192) + j) + 1)
------array access v.addr (((i * 8192) + j) + 1)
------array access v.addr ((((i + 1) * 8192) + j) + 1)
------array access v.addr ((((i + 2) * 8192) + j) + 1)
------array access q.addr (((i + 1) * 8192) + j)
------array access p.addr (((i + 1) * 8192) + j)
------array access q.addr ((((i + 1) * 8192) + j) + 1)
----array access u.addr ((((i + 1) * 8192) + 8192) - 1)
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p.addr (((((i + 1) * 8192) + 8192) - 2) - j)
------array access u.addr (((((i + 1) * 8192) + 8192) - 1) - j)
------array access q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
------array access u.addr (((((i + 1) * 8192) + 8192) - 2) - j)

Finish analysis loops */ 
/* # of Out-most Loops: 2 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 81
------Sample number: 6707
------Sample number: 6707
----Sample number: 81
------Sample number: 6707
------Sample number: 6707
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* v_addr0	98285 */
/* p_addr4	98285 */
/* u_addr6	98285 */
/* u_addr7	98285 */
/* p_addr1	98285 */
/* v_addr2	98285 */
/* q_addr3	98285 */
/* q_addr11	98285 */
/* v_addr12	98285 */
/* p_addr13	98285 */
/* p_addr5	98285 */
/* v_addr16	98285 */
/* u_addr17	98285 */
/* p_addr18	98285 */
/* u_addr19	98285 */
/* q_addr20	98285 */
/* p_addr21	98285 */
/* u_addr8	98285 */
/* q_addr9	98285 */
/* p_addr10	98285 */
/* q_addr28	98285 */
/* u_addr29	98285 */
/* p_addr30	98285 */
/* u_addr31	98285 */
/* q_addr32	98285 */
/* u_addr33	98285 */
/* v_addr14	98285 */
/* q_addr15	98285 */
/* p_addr22	98285 */
/* v_addr23	98285 */
/* v_addr24	98285 */
/* v_addr25	98285 */
/* q_addr26	98285 */
/* p_addr27	98285 */
#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <thread>
#include <vector>
#ifdef PAPI_TIMER
#  include "papi_timer.h"
#endif
using namespace std;
using namespace placeholders;
map<int, double> per_threadcnt_expected_rt;
map<int, map<uint64_t, double>> per_tcount_RT;
map<int, map<uint64_t, double>> privateRT;
map<int, double> per_tcount_MR;
map<uint64_t, double> RT;
map<uint64_t, double> MR;
uint64_t sample_sum = 0;
double total_private_stat = 0.0;
double total_reuse = 0.0;
double total_neighbor = 0.0;
double total_scale = 0.0;
double total_fold = 0.0;
double total_interchunk = 0.0;
double share_reuse = 0.0;
void privateRTHistoCal(int tid, map<int, map<uint64_t, double>> &rth, uint64_t rt, double val) {
    if (rth.find(tid) == rth.end()) { 
        rth[tid][rt] = val;
    } else {
        rth[tid][rt] += val;
    }
    return;
}
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
bool isCompleteChunk(uint64_t is, int thread_cnt) {
    return ((is % (CHUNK_SIZE * thread_cnt)) == 0);
}
int getChunkNum(uint64_t is, int thread_cnt) {
    if (is % (CHUNK_SIZE * thread_cnt) != 0) {
        return is / (CHUNK_SIZE * thread_cnt) + 1;
    }
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
    int sink_neighbor_delta = 0;
    if (!is_in_same_loop || getChunkID(i_src, thread_cnt) != getChunkID(i_sink, thread_cnt)) {
#ifdef DEBUG
        cout << "Inter Chunk Reuse" << endl;
        ;
        cout << "rt " << rt << endl;
#endif
        type = 0; // code for inter chunk reuse
        parallel_rt = rt * thread_cnt - CHUNK_SIZE * thread_cnt * (lsrc*(thread_cnt - tsrc) + lsink * tsink) + CHUNK_SIZE * thread_cnt * lsrc - (thread_cnt - 1) * middle_accesses + dT;
        if (dT != 0) {
            share_reuse += 1.0;
        } else {
            uint64_t private_rt = (CHUNK_SIZE - getThreadLocalPos(i_src)) * lsrc + (getThreadLocalPos(i_sink) + 1) * lsink + (middle_accesses / thread_cnt );
            privateRTHistoCal(tsrc, privateRT, private_rt, 1.0);
            total_private_stat += 1.0;
        }
        return parallel_rt;
    } else if (!is_normal_ref) {
        /* intra chunk reuse */
#ifdef DEBUG
        cout << "Neighboring Effect" << endl;
#endif
        if (tsrc == tsink) {
            privateRTHistoCal(tsrc, privateRT, rt, 1.0);
            total_private_stat += 1.0;
        }
        int tsrc_neighbor = search_src_candidate_neighbor(i_src, thread_cnt, srcAddrCal);
        int tsink_neighbor = search_sink_candidate_neighbor(i_sink, thread_cnt, sinkAddrCal);
        if (tsrc_neighbor >= 0) {
#ifdef DEBUG
            cout << "Find sink in src neighbor at " << tsrc_neighbor << endl;
            cout << "Neighbor Effect: " << tsrc_neighbor - tsrc << endl;
#endif
            total_neighbor += 1.0;
            type = 1; // code for src neighboring effect
            return tsrc_neighbor - tsrc;
        } else if (tsink_neighbor >= 0 && tsink_neighbor != tsink) {
#ifdef DEBUG
            cout << "Find sink in sink neighbor at " << tsink_neighbor << endl;
            cout << "Neighbor Effect: " << rt * thread_cnt + tsink_neighbor - tsink << endl;
#endif
            if (getChunkID(i_src, thread_cnt) == getChunkID(i_sink, thread_cnt)) {
                type = 2; // code for sink neighboring effect
                sink_neighbor_delta = tsink_neighbor - tsink;
                total_neighbor += 1.0;
                share_reuse += 1.0;
            }
        }
    }
    if (getChunkID(i_src, thread_cnt) == getChunkID(i_sink, thread_cnt)) {
        /* same thread -- scaling effect */
        if (dT == 0) {
#ifdef DEBUG
            cout << "Scaling Effect" << endl;
#endif
            total_scale += 1.0;
            parallel_rt = rt * thread_cnt;
            type = 3; // code for scaling effect
            privateRTHistoCal(tsrc, privateRT, rt, 1.0);
            total_private_stat += 1.0;
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
    return parallel_rt + sink_neighbor_delta;
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
double RTtoMR_AET_C(map<uint64_t, double> &rt, uint64_t cache_size) {
    map<uint64_t, double> P;
    double total_num_RT = 0;
    uint64_t max_RT = 0;
    for ( map<uint64_t, double>::reverse_iterator it = rt.rbegin(), eit = rt.rend(); it != eit; ++it) {
        total_num_RT += it->second;
        if (max_RT < it->first) {
            max_RT = it->first;
        }
    }
    double accumulate_num_RT = 0;
    for ( map<uint64_t, double>::reverse_iterator it = rt.rbegin(), eit = rt.rend(); it != eit; ++it) {
        P[it->first] = accumulate_num_RT / total_num_RT;
        accumulate_num_RT += it->second;
    }
    P[0] = 1;
    double sum_P = 0.0;
    uint64_t aet = 0, prev_aet = 0;
    while (sum_P < cache_size && aet <= max_RT) {
        if (P.find(aet) != P.end()) {
            sum_P += P[aet];
            prev_aet = aet;
        } else {
            sum_P += P[prev_aet];
        }
        aet++;
    }
    return P[prev_aet];
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
    cout << "Total Shared Reuses: " << share_reuse / total_reuse << endl;
    cout << "Total Reuses: " << total_reuse << endl;
    return;
}
/* Array v_addr	i */ 
/* i */
/* v_addr ((0 + i) + 1) 0 */
int calAddrv_addr0( int i) {
    int result = (((0 + i) + 1)) * 8 / 64;
    return result;
}
/* Array p_addr	i */ 
/* i */
/* p_addr (((i + 1) * 8192) + 0) 1 */
int calAddrp_addr1( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array v_addr	i */ 
/* i */
/* v_addr (0 + (i + 1)) 2 */
int calAddrv_addr2( int i) {
    int result = ((0 + (i + 1))) * 8 / 64;
    return result;
}
/* Array q_addr	i */ 
/* i */
/* q_addr (((i + 1) * 8192) + 0) 3 */
int calAddrq_addr3( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + j) 4 */
int calAddrp_addr4( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr ((((i + 1) * 8192) + j) + 1) 5 */
int calAddrp_addr5( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array u_addr	j i */ 
/* i */
/* i */
/* i, j */
/* u_addr (((j + 1) * 8192) + i) 6 */
int calAddru_addr6( int i, int j) {
    int result = ((((j + 1) * 8192) + i)) * 8 / 64;
    return result;
}
/* Array u_addr	j i */ 
/* i */
/* i */
/* i, j */
/* u_addr ((((j + 1) * 8192) + i) + 1) 7 */
int calAddru_addr7( int i, int j) {
    int result = (((((j + 1) * 8192) + i) + 1)) * 8 / 64;
    return result;
}
/* Array u_addr	j i */ 
/* i */
/* i */
/* i, j */
/* u_addr ((((j + 1) * 8192) + i) + 2) 8 */
int calAddru_addr8( int i, int j) {
    int result = (((((j + 1) * 8192) + i) + 2)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr (((i + 1) * 8192) + j) 9 */
int calAddrq_addr9( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + j) 10 */
int calAddrp_addr10( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr ((((i + 1) * 8192) + j) + 1) 11 */
int calAddrq_addr11( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array v_addr	i */ 
/* i */
/* v_addr (67100672 + (i + 1)) 12 */
int calAddrv_addr12( int i) {
    int result = ((67100672 + (i + 1))) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + (8190 - j)) 13 */
int calAddrp_addr13( int i, int j) {
    int result = ((((i + 1) * 8192) + (8190 - j))) * 8 / 64;
    return result;
}
/* Array v_addr	j i */ 
/* i */
/* i */
/* i, j */
/* v_addr ((((8191 - j) * 8192) + i) + 1) 14 */
int calAddrv_addr14( int i, int j) {
    int result = (((((8191 - j) * 8192) + i) + 1)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr (((((i + 1) * 8192) + 8192) - 2) - j) 15 */
int calAddrq_addr15( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
/* Array v_addr	j i */ 
/* i */
/* i */
/* i, j */
/* v_addr (((8190 - j) * 8192) + (i + 1)) 16 */
int calAddrv_addr16( int i, int j) {
    int result = ((((8190 - j) * 8192) + (i + 1))) * 8 / 64;
    return result;
}
/* Array u_addr	i */ 
/* i */
/* u_addr (((i + 1) * 8192) + 0) 17 */
int calAddru_addr17( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array p_addr	i */ 
/* i */
/* p_addr (((i + 1) * 8192) + 0) 18 */
int calAddrp_addr18( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array u_addr	i */ 
/* i */
/* u_addr (((i + 1) * 8192) + 0) 19 */
int calAddru_addr19( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array q_addr	i */ 
/* i */
/* q_addr (((i + 1) * 8192) + 0) 20 */
int calAddrq_addr20( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + j) 21 */
int calAddrp_addr21( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr ((((i + 1) * 8192) + j) + 1) 22 */
int calAddrp_addr22( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array v_addr	i j */ 
/* i */
/* v_addr (((i * 8192) + j) + 1) 23 */
int calAddrv_addr23( int i, int j) {
    int result = ((((i * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array v_addr	i j */ 
/* i */
/* v_addr ((((i + 1) * 8192) + j) + 1) 24 */
int calAddrv_addr24( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array v_addr	i j */ 
/* i */
/* v_addr ((((i + 2) * 8192) + j) + 1) 25 */
int calAddrv_addr25( int i, int j) {
    int result = (((((i + 2) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr (((i + 1) * 8192) + j) 26 */
int calAddrq_addr26( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + j) 27 */
int calAddrp_addr27( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr ((((i + 1) * 8192) + j) + 1) 28 */
int calAddrq_addr28( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array u_addr	i */ 
/* i */
/* u_addr ((((i + 1) * 8192) + 8192) - 1) 29 */
int calAddru_addr29( int i) {
    int result = (((((i + 1) * 8192) + 8192) - 1)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((((i + 1) * 8192) + 8192) - 2) - j) 30 */
int calAddrp_addr30( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
/* Array u_addr	i j */ 
/* i */
/* u_addr (((((i + 1) * 8192) + 8192) - 1) - j) 31 */
int calAddru_addr31( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 1) - j)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr (((((i + 1) * 8192) + 8192) - 2) - j) 32 */
int calAddrq_addr32( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
/* Array u_addr	i j */ 
/* i */
/* u_addr (((((i + 1) * 8192) + 8192) - 2) - j) 33 */
int calAddru_addr33( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
void ref_v_addr0(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr0( i) == calAddrv_addr0(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr0 --> v_addr0] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr0;" << "(" << i_Start<< ");" << "v_addr0;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr2( i) == calAddrv_addr0(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr0 --> v_addr2] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr0;" << "(" << i_Start<< ");" << "v_addr2;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr12( i) == calAddrv_addr0(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr0 --> v_addr12] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr0;" << "(" << i_Start<< ");" << "v_addr12;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr14( i, j) == calAddrv_addr0(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr0 --> v_addr14] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr0;" << "(" << i_Start<< ");" << "v_addr14;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr16( i, j) == calAddrv_addr0(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr0 --> v_addr16] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr0;" << "(" << i_Start<< ");" << "v_addr16;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr23( i, j) == calAddrv_addr0(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr0 --> v_addr23] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr0;" << "(" << i_Start<< ");" << "v_addr23;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr24( i, j) == calAddrv_addr0(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr0 --> v_addr24] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr0;" << "(" << i_Start<< ");" << "v_addr24;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr25( i, j) == calAddrv_addr0(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr0 --> v_addr25] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr0;" << "(" << i_Start<< ");" << "v_addr25;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr4(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr1( i) == calAddrp_addr4(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr4 --> p_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr1;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr4( i, j) == calAddrp_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr4;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( i, j) == calAddrp_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr10( i, j) == calAddrp_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr10;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr13( i, j) == calAddrp_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr13;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr4(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr4 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr6(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr6( i, j) == calAddru_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr6 --> u_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr6;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr7( i, j) == calAddru_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr6 --> u_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr7;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr8( i, j) == calAddru_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr6 --> u_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr8;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr17( i) == calAddru_addr6(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr6 --> u_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr17;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr19( i) == calAddru_addr6(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr6 --> u_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr19;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr29( i) == calAddru_addr6(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr6 --> u_addr29] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr29;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr31( i, j) == calAddru_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr6 --> u_addr31] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr31;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr33( i, j) == calAddru_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr6 --> u_addr33] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr33;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
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
void ref_u_addr7(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr6( i, j) == calAddru_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr7 --> u_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr6;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr7( i, j) == calAddru_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr7 --> u_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr7;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr8( i, j) == calAddru_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr7 --> u_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr8;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr17( i) == calAddru_addr7(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr7 --> u_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr17;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr19( i) == calAddru_addr7(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr7 --> u_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr19;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr29( i) == calAddru_addr7(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr7 --> u_addr29] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr29;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr31( i, j) == calAddru_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr7 --> u_addr31] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr31;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr33( i, j) == calAddru_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr7 --> u_addr33] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr33;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
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
void ref_p_addr1(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr1( i) == calAddrp_addr1(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr1 --> p_addr1] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr1;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr4( i, j) == calAddrp_addr1(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr4] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr4;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( i, j) == calAddrp_addr1(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr5] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr10( i, j) == calAddrp_addr1(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr10] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr10;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr13( i, j) == calAddrp_addr1(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr13] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr13;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr1(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr1 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr1(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr1(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr1(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr1(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr1;" << "(" << i_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr2(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr0( i) == calAddrv_addr2(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr2 --> v_addr0] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr2;" << "(" << i_Start<< ");" << "v_addr0;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr2( i) == calAddrv_addr2(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr2 --> v_addr2] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr2;" << "(" << i_Start<< ");" << "v_addr2;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr12( i) == calAddrv_addr2(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr2 --> v_addr12] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr2;" << "(" << i_Start<< ");" << "v_addr12;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr14( i, j) == calAddrv_addr2(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr2 --> v_addr14] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr2;" << "(" << i_Start<< ");" << "v_addr14;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr16( i, j) == calAddrv_addr2(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr2 --> v_addr16] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr2;" << "(" << i_Start<< ");" << "v_addr16;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr23( i, j) == calAddrv_addr2(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr2 --> v_addr23] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr2;" << "(" << i_Start<< ");" << "v_addr23;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr24( i, j) == calAddrv_addr2(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr2 --> v_addr24] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr2;" << "(" << i_Start<< ");" << "v_addr24;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr25( i, j) == calAddrv_addr2(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr2 --> v_addr25] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr2;" << "(" << i_Start<< ");" << "v_addr25;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr3(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr3( i) == calAddrq_addr3(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr3 --> q_addr3] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr3;" << "(" << i_Start<< ");" << "q_addr3;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr9( i, j) == calAddrq_addr3(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr3 --> q_addr9] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr3;" << "(" << i_Start<< ");" << "q_addr9;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr11( i, j) == calAddrq_addr3(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr3 --> q_addr11] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr3;" << "(" << i_Start<< ");" << "q_addr11;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr15( i, j) == calAddrq_addr3(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr3 --> q_addr15] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr3;" << "(" << i_Start<< ");" << "q_addr15;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr20( i) == calAddrq_addr3(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr3 --> q_addr20] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr3;" << "(" << i_Start<< ");" << "q_addr20;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr26( i, j) == calAddrq_addr3(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr3 --> q_addr26] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr3;" << "(" << i_Start<< ");" << "q_addr26;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr28( i, j) == calAddrq_addr3(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr3 --> q_addr28] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr3;" << "(" << i_Start<< ");" << "q_addr28;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr32( i, j) == calAddrq_addr3(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr3 --> q_addr32] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr3;" << "(" << i_Start<< ");" << "q_addr32;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr11(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr3( i) == calAddrq_addr11(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr11 --> q_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr11;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr3;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr9( i, j) == calAddrq_addr11(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr11 --> q_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr11;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr9;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr11( i, j) == calAddrq_addr11(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr11 --> q_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr11;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr11;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr15( i, j) == calAddrq_addr11(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr11 --> q_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr11;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr15;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr20( i) == calAddrq_addr11(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr11 --> q_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr11;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr20;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr26( i, j) == calAddrq_addr11(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr11 --> q_addr26] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr11;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr26;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr28( i, j) == calAddrq_addr11(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr11 --> q_addr28] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr11;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr28;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr32( i, j) == calAddrq_addr11(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr11 --> q_addr32] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr11;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr32;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr12(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr0( i) == calAddrv_addr12(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr12 --> v_addr0] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr12;" << "(" << i_Start<< ");" << "v_addr0;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr2( i) == calAddrv_addr12(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr12 --> v_addr2] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr12;" << "(" << i_Start<< ");" << "v_addr2;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr12( i) == calAddrv_addr12(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr12 --> v_addr12] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr12;" << "(" << i_Start<< ");" << "v_addr12;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr14( i, j) == calAddrv_addr12(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr12 --> v_addr14] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr12;" << "(" << i_Start<< ");" << "v_addr14;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr16( i, j) == calAddrv_addr12(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr12 --> v_addr16] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr12;" << "(" << i_Start<< ");" << "v_addr16;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr23( i, j) == calAddrv_addr12(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr12 --> v_addr23] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr12;" << "(" << i_Start<< ");" << "v_addr23;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr24( i, j) == calAddrv_addr12(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr12 --> v_addr24] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr12;" << "(" << i_Start<< ");" << "v_addr24;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr25( i, j) == calAddrv_addr12(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr12 --> v_addr25] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr12;" << "(" << i_Start<< ");" << "v_addr25;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr13(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr1( i) == calAddrp_addr13(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr13 --> p_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr1;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr4( i, j) == calAddrp_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr4;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( i, j) == calAddrp_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr10( i, j) == calAddrp_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr10;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr13( i, j) == calAddrp_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr13;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr13(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr13 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr5(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr1( i) == calAddrp_addr5(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr5 --> p_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr1;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr4( i, j) == calAddrp_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr4;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( i, j) == calAddrp_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr10( i, j) == calAddrp_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr10;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr13( i, j) == calAddrp_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr13;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr5(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr5 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr16(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr0( i) == calAddrv_addr16(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr16 --> v_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr0;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr2( i) == calAddrv_addr16(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr16 --> v_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr2;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr12( i) == calAddrv_addr16(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr16 --> v_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr12;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr14( i, j) == calAddrv_addr16(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr16 --> v_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr14;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr16( i, j) == calAddrv_addr16(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr16 --> v_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr16;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr23( i, j) == calAddrv_addr16(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr16 --> v_addr23] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr23;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr24( i, j) == calAddrv_addr16(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr16 --> v_addr24] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr24;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr25( i, j) == calAddrv_addr16(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr16 --> v_addr25] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr25;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr17(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr17( i) == calAddru_addr17(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr17 --> u_addr17] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr17;" << "(" << i_Start<< ");" << "u_addr17;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr19( i) == calAddru_addr17(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr17 --> u_addr19] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr17;" << "(" << i_Start<< ");" << "u_addr19;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr29( i) == calAddru_addr17(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr17 --> u_addr29] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr17;" << "(" << i_Start<< ");" << "u_addr29;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr31( i, j) == calAddru_addr17(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr17 --> u_addr31] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr17;" << "(" << i_Start<< ");" << "u_addr31;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr33( i, j) == calAddru_addr17(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr17 --> u_addr33] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr17;" << "(" << i_Start<< ");" << "u_addr33;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
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
void ref_p_addr18(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr18(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr18 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr18;" << "(" << i_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr18(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr18 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr18;" << "(" << i_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr18(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr18 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr18;" << "(" << i_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr18(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr18 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr18;" << "(" << i_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr18(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr18 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr18;" << "(" << i_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr19(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr17( i) == calAddru_addr19(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr19 --> u_addr17] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr19;" << "(" << i_Start<< ");" << "u_addr17;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr19( i) == calAddru_addr19(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr19 --> u_addr19] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr19;" << "(" << i_Start<< ");" << "u_addr19;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr29( i) == calAddru_addr19(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr19 --> u_addr29] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr19;" << "(" << i_Start<< ");" << "u_addr29;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr31( i, j) == calAddru_addr19(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr19 --> u_addr31] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr19;" << "(" << i_Start<< ");" << "u_addr31;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr33( i, j) == calAddru_addr19(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr19 --> u_addr33] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr19;" << "(" << i_Start<< ");" << "u_addr33;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
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
void ref_q_addr20(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr20( i) == calAddrq_addr20(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr20 --> q_addr20] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr20;" << "(" << i_Start<< ");" << "q_addr20;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr26( i, j) == calAddrq_addr20(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr20 --> q_addr26] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr20;" << "(" << i_Start<< ");" << "q_addr26;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr28( i, j) == calAddrq_addr20(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr20 --> q_addr28] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr20;" << "(" << i_Start<< ");" << "q_addr28;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr32( i, j) == calAddrq_addr20(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr20 --> q_addr32] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr20;" << "(" << i_Start<< ");" << "q_addr32;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr21(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr21(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr21 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr21;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr21(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr21 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr21;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr21(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr21 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr21;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr21(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr21 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr21;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr21(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr21 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr21;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr8(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr6( i, j) == calAddru_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr8 --> u_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr6;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr7( i, j) == calAddru_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr8 --> u_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr7;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr8( i, j) == calAddru_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr8 --> u_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr8;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr17( i) == calAddru_addr8(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr8 --> u_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr17;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr19( i) == calAddru_addr8(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr8 --> u_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr19;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr29( i) == calAddru_addr8(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr8 --> u_addr29] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr29;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr31( i, j) == calAddru_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr8 --> u_addr31] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr31;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr33( i, j) == calAddru_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr8 --> u_addr33] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr33;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
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
void ref_q_addr9(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr3( i) == calAddrq_addr9(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr9 --> q_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr3;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr9( i, j) == calAddrq_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr9 --> q_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr9;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr11( i, j) == calAddrq_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr9 --> q_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr11;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr15( i, j) == calAddrq_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr9 --> q_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr15;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr20( i) == calAddrq_addr9(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr9 --> q_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr20;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr26( i, j) == calAddrq_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr9 --> q_addr26] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr26;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr28( i, j) == calAddrq_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr9 --> q_addr28] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr28;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr32( i, j) == calAddrq_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr9 --> q_addr32] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr32;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr10(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr1( i) == calAddrp_addr10(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr10 --> p_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr1;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr4( i, j) == calAddrp_addr10(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr4;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( i, j) == calAddrp_addr10(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr10( i, j) == calAddrp_addr10(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr10;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr13( i, j) == calAddrp_addr10(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr13;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr10(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr10 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr10(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr10(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr10(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr10(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr28(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr20( i) == calAddrq_addr28(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr28 --> q_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr28;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr20;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr26( i, j) == calAddrq_addr28(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr28 --> q_addr26] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr28;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr26;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr28( i, j) == calAddrq_addr28(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr28 --> q_addr28] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr28;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr28;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr32( i, j) == calAddrq_addr28(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr28 --> q_addr32] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr28;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr32;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr29(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 81;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr17( i) == calAddru_addr29(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr29 --> u_addr17] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr29;" << "(" << i_Start<< ");" << "u_addr17;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr19( i) == calAddru_addr29(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr29 --> u_addr19] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr29;" << "(" << i_Start<< ");" << "u_addr19;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr29( i) == calAddru_addr29(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr29 --> u_addr29] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr29;" << "(" << i_Start<< ");" << "u_addr29;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr31( i, j) == calAddru_addr29(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr29 --> u_addr31] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr29;" << "(" << i_Start<< ");" << "u_addr31;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr33( i, j) == calAddru_addr29(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr29 --> u_addr33] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr29;" << "(" << i_Start<< ");" << "u_addr33;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
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
void ref_p_addr30(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr30(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr30 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr30;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr30(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr30 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr30;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr30(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr30 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr30;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr30(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr30 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr30;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr30(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr30 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr30;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr31(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr17( i) == calAddru_addr31(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr31 --> u_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr31;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr17;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr19( i) == calAddru_addr31(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr31 --> u_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr31;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr19;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr29( i) == calAddru_addr31(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr31 --> u_addr29] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr31;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr29;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr31( i, j) == calAddru_addr31(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr31 --> u_addr31] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr31;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr31;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr33( i, j) == calAddru_addr31(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr31 --> u_addr33] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr31;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr33;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
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
void ref_q_addr32(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr20( i) == calAddrq_addr32(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr32 --> q_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr32;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr20;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr26( i, j) == calAddrq_addr32(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr32 --> q_addr26] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr32;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr26;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr28( i, j) == calAddrq_addr32(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr32 --> q_addr28] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr32;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr28;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr32( i, j) == calAddrq_addr32(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr32 --> q_addr32] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr32;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr32;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr33(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr17( i) == calAddru_addr33(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr33 --> u_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr33;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr17;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr19( i) == calAddru_addr33(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr33 --> u_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr33;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr19;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddru_addr29( i) == calAddru_addr33(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[u_addr33 --> u_addr29] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "u_addr33;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr29;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr31( i, j) == calAddru_addr33(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr33 --> u_addr31] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr33;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr31;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr33( i, j) == calAddru_addr33(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u_addr33 --> u_addr33] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u_addr33;" << "(" << i_Start<< ", " << j_Start<< ");" << "u_addr33;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr14(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr0( i) == calAddrv_addr14(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr14 --> v_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr14;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr0;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr2( i) == calAddrv_addr14(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr14 --> v_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr14;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr2;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrv_addr12( i) == calAddrv_addr14(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[v_addr14 --> v_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "v_addr14;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr12;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr14( i, j) == calAddrv_addr14(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr14 --> v_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr14;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr14;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr16( i, j) == calAddrv_addr14(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr14 --> v_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr14;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr16;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr23( i, j) == calAddrv_addr14(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr14 --> v_addr23] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr14;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr23;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr24( i, j) == calAddrv_addr14(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr14 --> v_addr24] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr14;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr24;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr25( i, j) == calAddrv_addr14(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr14 --> v_addr25] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr14;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr25;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr15(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr3( i) == calAddrq_addr15(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr15 --> q_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr15;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr3;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr9( i, j) == calAddrq_addr15(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr15 --> q_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr15;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr9;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr11( i, j) == calAddrq_addr15(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr15 --> q_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr15;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr11;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr15( i, j) == calAddrq_addr15(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr15 --> q_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr15;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr15;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr20( i) == calAddrq_addr15(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr15 --> q_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr15;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr20;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr26( i, j) == calAddrq_addr15(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr15 --> q_addr26] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr15;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr26;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr28( i, j) == calAddrq_addr15(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr15 --> q_addr28] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr15;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr28;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr32( i, j) == calAddrq_addr15(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8190 % (CHUNK_SIZE * thread_cnt))) * 98285 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr15 --> q_addr32] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr15;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr32;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr22(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr22(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr22 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr22;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr22(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr22 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr22;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr22(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr22 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr22;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr22(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr22 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr22;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr22(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr22 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr22;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr23(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr23( i, j) == calAddrv_addr23(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr23 --> v_addr23] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr23;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr23;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr24( i, j) == calAddrv_addr23(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr23 --> v_addr24] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr23;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr24;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr25( i, j) == calAddrv_addr23(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr23 --> v_addr25] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr23;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr25;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr24(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr23( i, j) == calAddrv_addr24(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr24 --> v_addr23] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr24;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr23;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr24( i, j) == calAddrv_addr24(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr24 --> v_addr24] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr24;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr24;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr25( i, j) == calAddrv_addr24(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr24 --> v_addr25] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr24;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr25;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr25(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr23( i, j) == calAddrv_addr25(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr25 --> v_addr23] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr25;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr23;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr24( i, j) == calAddrv_addr25(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr25 --> v_addr24] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr25;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr24;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr25( i, j) == calAddrv_addr25(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v_addr25 --> v_addr25] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v_addr25;" << "(" << i_Start<< ", " << j_Start<< ");" << "v_addr25;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr26(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr20( i) == calAddrq_addr26(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[q_addr26 --> q_addr20] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "q_addr26;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr20;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr26( i, j) == calAddrq_addr26(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr26 --> q_addr26] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr26;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr26;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr28( i, j) == calAddrq_addr26(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr26 --> q_addr28] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr26;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr28;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr32( i, j) == calAddrq_addr26(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[q_addr26 --> q_addr32] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "q_addr26;" << "(" << i_Start<< ", " << j_Start<< ");" << "q_addr32;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr27(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6707;) {
SAMPLE:
        int i_Start = rand() % (8190 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8190 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrp_addr18( i) == calAddrp_addr27(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[p_addr27 --> p_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "p_addr27;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr18;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr21( i, j) == calAddrp_addr27(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr27 --> p_addr21] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr27;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr21;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr22( i, j) == calAddrp_addr27(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr27 --> p_addr22] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr27;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr22;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr27( i, j) == calAddrp_addr27(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr27 --> p_addr27] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr27;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr27;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr30( i, j) == calAddrp_addr27(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 98285, 98285, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[p_addr27 --> p_addr30] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "p_addr27;" << "(" << i_Start<< ", " << j_Start<< ");" << "p_addr30;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void mergePrivateRT(map<int, map<uint64_t, double>> prt) {
    for(map<int, map<uint64_t, double>>::iterator it = prt.begin(); it != prt.end(); ++it) {
        map<uint64_t, double> tmp = it->second;
        for(map<uint64_t, double>::iterator iit = tmp.begin(); iit != tmp.end(); ++iit) {
            rtHistoCal(RT, iit->first, iit->second);
        }
    }
    return;
}
void generate_per_thread_reuse(int thread_cnt) {
    map<uint64_t, double> tmpRT;
    double accumulate_num_RT = 0.0;
    double expected_rt = 0.0;
    ref_v_addr0(thread_cnt, tmpRT);
    ref_p_addr4(thread_cnt, tmpRT);
    ref_u_addr6(thread_cnt, tmpRT);
    ref_u_addr7(thread_cnt, tmpRT);
    ref_p_addr1(thread_cnt, tmpRT);
    ref_v_addr2(thread_cnt, tmpRT);
    ref_q_addr3(thread_cnt, tmpRT);
    ref_q_addr11(thread_cnt, tmpRT);
    ref_v_addr12(thread_cnt, tmpRT);
    ref_p_addr13(thread_cnt, tmpRT);
    ref_p_addr5(thread_cnt, tmpRT);
    ref_v_addr16(thread_cnt, tmpRT);
    ref_u_addr17(thread_cnt, tmpRT);
    ref_p_addr18(thread_cnt, tmpRT);
    ref_u_addr19(thread_cnt, tmpRT);
    ref_q_addr20(thread_cnt, tmpRT);
    ref_p_addr21(thread_cnt, tmpRT);
    ref_u_addr8(thread_cnt, tmpRT);
    ref_q_addr9(thread_cnt, tmpRT);
    ref_p_addr10(thread_cnt, tmpRT);
    ref_q_addr28(thread_cnt, tmpRT);
    ref_u_addr29(thread_cnt, tmpRT);
    ref_p_addr30(thread_cnt, tmpRT);
    ref_u_addr31(thread_cnt, tmpRT);
    ref_q_addr32(thread_cnt, tmpRT);
    ref_u_addr33(thread_cnt, tmpRT);
    ref_v_addr14(thread_cnt, tmpRT);
    ref_q_addr15(thread_cnt, tmpRT);
    ref_p_addr22(thread_cnt, tmpRT);
    ref_v_addr23(thread_cnt, tmpRT);
    ref_v_addr24(thread_cnt, tmpRT);
    ref_v_addr25(thread_cnt, tmpRT);
    ref_q_addr26(thread_cnt, tmpRT);
    ref_p_addr27(thread_cnt, tmpRT);
    per_tcount_RT[thread_cnt] = tmpRT;
    /* iterate each tid and dump its LLC miss ratio (13.75MB) */
    per_tcount_MR[thread_cnt] = RTtoMR_AET_C(tmpRT, 2285280);
}
/* Compute expected reuse for current thread number */
void compute_expected_reuse(int thread_cnt) {
    double accumulate_num_RT = 0.0;
    double expected_rt = 0.0;
    map<uint64_t, double> tmpRT = per_tcount_RT[thread_cnt];
    for ( map<uint64_t, double>::iterator it = tmpRT.begin(), eit = tmpRT.end(); it != eit; ++it) {
        accumulate_num_RT += it->second;
    }
    for ( map<uint64_t, double>::iterator it = tmpRT.begin(), eit = tmpRT.end(); it != eit; ++it) {
        expected_rt += (it->first * it->second / accumulate_num_RT);
    }
#ifdef DEBUG
    cout << thread_cnt << ", " << expected_rt << endl;
#endif
    per_threadcnt_expected_rt[thread_cnt] = expected_rt;
}
int main() {
    int tlb = THREAD_LB;
    int tub = THREAD_UB;
    /* metadata to derive best thread */
    double min_expected_rt = 0.0;
    int min_tnum = tlb;
    vector<thread> thread_vec;
    for (int t = tlb; t <= tub; t++) {
        /* Currently we consider all even number threads only */
        if ( ceil(log2(t)) != floor(log2(t)) ) { continue; }
        generate_per_thread_reuse(t);
        /* iterate each tid and dump its private L2 miss ratio (1MB) */
        for (int id = 0; id < t; id++) {
            cout << RTtoMR_AET_C(privateRT[id], 16384) << " | ";
        }
        cout << endl;        thread_vec.push_back(thread(compute_expected_reuse, t));
    }
    /* compute expected rt for each thread count */
    for (vector<thread>::iterator threadit = thread_vec.begin(); threadit != thread_vec.end(); ++threadit) {
        (*threadit).join();
    }
    for (map<int, double>::iterator mit = per_threadcnt_expected_rt.begin(); mit != per_threadcnt_expected_rt.end(); ++mit) {
        cout << mit->first << ", " << mit->second << endl;
        if (mit == per_threadcnt_expected_rt.begin()) {
;            min_expected_rt = mit->second;
            continue;
        }
        if (min_expected_rt > mit->second) {
            min_tnum = mit->first;
            min_expected_rt = mit->second;
        }
    }
    cout << min_tnum << endl;
    return 0;
}
 /* Analyze function: adi */ 
