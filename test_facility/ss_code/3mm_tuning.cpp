
 /* Start to analysis array index
Array index info: Total number of references: 15
E.addr ((i * 1024) + j)
E.addr ((i * 1024) + j)
E.addr ((i * 1024) + j)
C.addr ((i * 1024) + k)
D.addr ((k * 1024) + j)
A.addr ((i * 1024) + k)
B.addr ((k * 1024) + j)
G.addr ((i * 1024) + j)
E.addr ((i * 1024) + k)
F.addr ((k * 1024) + j)
F.addr ((i * 1024) + j)
F.addr ((i * 1024) + j)
F.addr ((i * 1024) + j)
G.addr ((i * 1024) + j)
G.addr ((i * 1024) + j)
BC Array cost info: Total number of arrays: 7
A.addr 2
B.addr 2
C.addr 2
D.addr 2
E.addr 8
F.addr 8
G.addr 6
BC Average cost: 4.285714e+00

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %ni
i32 %nj
i32 %nk
i32 %nl
i32 %nm
double* %E
double* %A
double* %B
double* %F
double* %C
double* %D
double* %G

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access E.addr ((i * 1024) + j)
------k
------Loop Bound: (0, 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 1024) + k)
--------array access B.addr ((k * 1024) + j)
--------array access E.addr ((i * 1024) + j)
--------array access E.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access F.addr ((i * 1024) + j)
------k
------Loop Bound: (0, 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access C.addr ((i * 1024) + k)
--------array access D.addr ((k * 1024) + j)
--------array access F.addr ((i * 1024) + j)
--------array access F.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access G.addr ((i * 1024) + j)
------k
------Loop Bound: (0, 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access E.addr ((i * 1024) + k)
--------array access F.addr ((k * 1024) + j)
--------array access G.addr ((i * 1024) + j)
--------array access G.addr ((i * 1024) + j)

Finish analysis loops */ 
/* # of Out-most Loops: 3 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 10
------Sample number: 104
--------Sample number: 1073
----Sample number: 10
------Sample number: 104
--------Sample number: 1073
----Sample number: 10
------Sample number: 104
--------Sample number: 1073
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* E_addr3	4195328 */
/* E_addr4	4195328 */
/* E_addr0	4195328 */
/* C_addr6	4195328 */
/* D_addr7	4195328 */
/* A_addr1	4195328 */
/* B_addr2	4195328 */
/* G_addr10	4195328 */
/* E_addr11	4195328 */
/* F_addr12	4195328 */
/* F_addr5	4195328 */
/* F_addr8	4195328 */
/* F_addr9	4195328 */
/* G_addr13	4195328 */
/* G_addr14	4195328 */
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
/* Array E_addr	i j */ 
/* i */
/* E_addr ((i * 1024) + j) 0 */
int calAddrE_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array A_addr	i k */ 
/* i */
/* A_addr ((i * 1024) + k) 1 */
int calAddrA_addr1( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* Array B_addr	k j */ 
/* i */
/* i */
/* B_addr ((k * 1024) + j) 2 */
int calAddrB_addr2( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
/* Array E_addr	i j */ 
/* i */
/* E_addr ((i * 1024) + j) 3 */
int calAddrE_addr3( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array E_addr	i j */ 
/* i */
/* E_addr ((i * 1024) + j) 4 */
int calAddrE_addr4( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array F_addr	i j */ 
/* i */
/* F_addr ((i * 1024) + j) 5 */
int calAddrF_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array C_addr	i k */ 
/* i */
/* C_addr ((i * 1024) + k) 6 */
int calAddrC_addr6( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* Array D_addr	k j */ 
/* i */
/* i */
/* D_addr ((k * 1024) + j) 7 */
int calAddrD_addr7( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
/* Array F_addr	i j */ 
/* i */
/* F_addr ((i * 1024) + j) 8 */
int calAddrF_addr8( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array F_addr	i j */ 
/* i */
/* F_addr ((i * 1024) + j) 9 */
int calAddrF_addr9( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array G_addr	i j */ 
/* i */
/* G_addr ((i * 1024) + j) 10 */
int calAddrG_addr10( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array E_addr	i k */ 
/* i */
/* E_addr ((i * 1024) + k) 11 */
int calAddrE_addr11( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* Array F_addr	k j */ 
/* i */
/* i */
/* F_addr ((k * 1024) + j) 12 */
int calAddrF_addr12( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
/* Array G_addr	i j */ 
/* i */
/* G_addr ((i * 1024) + j) 13 */
int calAddrG_addr13( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array G_addr	i j */ 
/* i */
/* G_addr ((i * 1024) + j) 14 */
int calAddrG_addr14( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_E_addr3(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrE_addr0( i, j) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr3, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr0, _1, j);
                            /* i 4195328 */
                            /* j 4097 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[E_addr3 --> E_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "E_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "E_addr0;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr3, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr3 --> E_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "E_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr4( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr4, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr3 --> E_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "E_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr11, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                if (isCompleteChunk(1024, thread_cnt)) {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                } else {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (1024 % (CHUNK_SIZE * thread_cnt))) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                }
                                middle_accesses += 4296015872;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr3 --> E_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr3;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "E_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        }
        }
EndSample:
        s++;
        }
}
void ref_E_addr4(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrE_addr0( i, j) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr4, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr0, _1, j);
                            /* i 4195328 */
                            /* j 4097 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[E_addr4 --> E_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "E_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "E_addr0;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr3, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr4 --> E_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "E_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr4( i, j, k) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr4, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr4 --> E_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "E_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr11, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                if (isCompleteChunk(1024, thread_cnt)) {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                } else {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (1024 % (CHUNK_SIZE * thread_cnt))) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                }
                                middle_accesses += 4296015872;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr4 --> E_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr4;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "E_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        }
        }
EndSample:
        s++;
        }
}
void ref_E_addr0(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrE_addr0( i, j) == calAddrE_addr0(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr0, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr0, _1, j);
                            /* i 4195328 */
                            /* j 4097 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[E_addr0 --> E_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "E_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "E_addr0;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr0, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr3, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr0 --> E_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "E_addr3;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr4( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr0, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr4, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr0 --> E_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "E_addr4;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr0, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr11, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                if (isCompleteChunk(1024, thread_cnt)) {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                } else {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (1024 % (CHUNK_SIZE * thread_cnt))) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                }
                                middle_accesses += 4296015872;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr0 --> E_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "E_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        }
        }
EndSample:
        s++;
        }
}
void ref_C_addr6(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrC_addr6( i, j, k) == calAddrC_addr6(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrC_addr6, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrC_addr6, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[C_addr6 --> C_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "C_addr6;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "C_addr6;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        }
        }
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
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
void ref_D_addr7(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr7( i, j, k) == calAddrD_addr7(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr7, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[D_addr7 --> D_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "D_addr7;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "D_addr7;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
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
void ref_A_addr1(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k=(k + 1)) {
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
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
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
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
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
void ref_B_addr2(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( i, j, k) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr2, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[B_addr2 --> B_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "B_addr2;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "B_addr2;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
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
void ref_G_addr10(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrG_addr10( i, j) == calAddrG_addr10(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrG_addr10, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrG_addr10, _1, j);
                            /* i 4195328 */
                            /* j 4097 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[G_addr10 --> G_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "G_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "G_addr10;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr13( i, j, k) == calAddrG_addr10(i_Start, j_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrG_addr10, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrG_addr13, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[G_addr10 --> G_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "G_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "G_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr14( i, j, k) == calAddrG_addr10(i_Start, j_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrG_addr10, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrG_addr14, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[G_addr10 --> G_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "G_addr10;" << "(" << i_Start<< ", " << j_Start<< ");" << "G_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_E_addr11(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr11(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrE_addr11, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrE_addr11, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[E_addr11 --> E_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "E_addr11;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "E_addr11;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        }
        }
EndSample:
        s++;
        }
}
void ref_F_addr12(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr12(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr12, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr12, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr12 --> F_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr12;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "F_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
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
void ref_F_addr5(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrF_addr5( i, j) == calAddrF_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr5, _1, j);
                            /* i 4195328 */
                            /* j 4097 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[F_addr5 --> F_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "F_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "F_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr8( i, j, k) == calAddrF_addr5(i_Start, j_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr5, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr8, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr5 --> F_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "F_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr9( i, j, k) == calAddrF_addr5(i_Start, j_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr5, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr9, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr5 --> F_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "F_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr5(i_Start, j_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr5, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr12, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                if (isCompleteChunk(1024, thread_cnt)) {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                } else {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (1024 % (CHUNK_SIZE * thread_cnt))) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                }
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr5 --> F_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "F_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
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
void ref_F_addr8(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrF_addr5( i, j) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr8, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr5, _1, j);
                            /* i 4195328 */
                            /* j 4097 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[F_addr8 --> F_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "F_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "F_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr8( i, j, k) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr8, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr8 --> F_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "F_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr9( i, j, k) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr9, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr8 --> F_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "F_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr12, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                if (isCompleteChunk(1024, thread_cnt)) {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                } else {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (1024 % (CHUNK_SIZE * thread_cnt))) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                }
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr8 --> F_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr8;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "F_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
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
void ref_F_addr9(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrF_addr5( i, j) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr9, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr5, _1, j);
                            /* i 4195328 */
                            /* j 4097 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[F_addr9 --> F_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "F_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "F_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr8( i, j, k) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr8, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr9 --> F_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "F_addr8;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr9( i, j, k) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr9, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr9 --> F_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "F_addr9;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int iLB6 = 0;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: false */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrF_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrF_addr12, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                if (isCompleteChunk(1024, thread_cnt)) {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                } else {
                                                    middle_accesses += ((getChunkNum(1024, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (1024 % (CHUNK_SIZE * thread_cnt))) * 4195328 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 4195328;
                                }
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[F_addr9 --> F_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "F_addr9;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "F_addr12;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
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
void ref_G_addr13(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrG_addr10( i, j) == calAddrG_addr13(i_Start, j_Start, k_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrG_addr13, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrG_addr10, _1, j);
                            /* i 4195328 */
                            /* j 4097 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[G_addr13 --> G_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "G_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "G_addr10;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                {
                int kLB8 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr13( i, j, k) == calAddrG_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrG_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrG_addr13, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[G_addr13 --> G_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "G_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "G_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr14( i, j, k) == calAddrG_addr13(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrG_addr13, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrG_addr14, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[G_addr13 --> G_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "G_addr13;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "G_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_G_addr14(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 1024) { goto SAMPLE; }
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        if (k_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 1024; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrG_addr10( i, j) == calAddrG_addr14(i_Start, j_Start, k_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrG_addr14, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrG_addr10, _1, j);
                            /* i 4195328 */
                            /* j 4097 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[G_addr14 --> G_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "G_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "G_addr10;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                {
                int kLB8 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr13( i, j, k) == calAddrG_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrG_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrG_addr13, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[G_addr14 --> G_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "G_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "G_addr13;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr14( i, j, k) == calAddrG_addr14(i_Start, j_Start, k_Start)) {
                                /* is_src_loop_outermost: 1 */
                                /* is_sink_loop_outermost: 1 */
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrG_addr14, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrG_addr14, _1, j, k);
                                /* i 4195328 */
                                /* j 4097 */
                                /* k 4 */
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0;
                                middle_accesses += 0;
#ifdef DEBUG
                                cout << " middle_access is " << middle_accesses << endl;
#endif
                                int reuse_type = -1;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                                rtHistoCal(RT, parallel_rt, 1.0);
#else
                                rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                    cout << "[G_addr14 --> G_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                                    // cout << "G_addr14;" << "(" << i_Start<< ", " << j_Start<< ", " << k_Start<< ");" << "G_addr14;" << "(" << i<< ", " << j<< ", " << k<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
    ref_E_addr3(thread_cnt, tmpRT);
    ref_E_addr4(thread_cnt, tmpRT);
    ref_E_addr0(thread_cnt, tmpRT);
    ref_C_addr6(thread_cnt, tmpRT);
    ref_D_addr7(thread_cnt, tmpRT);
    ref_A_addr1(thread_cnt, tmpRT);
    ref_B_addr2(thread_cnt, tmpRT);
    ref_G_addr10(thread_cnt, tmpRT);
    ref_E_addr11(thread_cnt, tmpRT);
    ref_F_addr12(thread_cnt, tmpRT);
    ref_F_addr5(thread_cnt, tmpRT);
    ref_F_addr8(thread_cnt, tmpRT);
    ref_F_addr9(thread_cnt, tmpRT);
    ref_G_addr13(thread_cnt, tmpRT);
    ref_G_addr14(thread_cnt, tmpRT);
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
        if ( t % 2 != 0 ) { continue; }
        generate_per_thread_reuse(t);
        /* iterate each tid and dump its private L2 miss ratio (1MB) */
        thread_vec.push_back(thread(compute_expected_reuse, t));
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
 /* Analyze function: mm2 */ 
