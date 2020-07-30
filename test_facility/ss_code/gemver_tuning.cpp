
 /* Start to analysis array index
Array index info: Total number of references: 17
A.addr ((i * 4096) + j)
A.addr ((i * 4096) + j)
u1.addr i
y.addr j
x.addr i
v1.addr j
u2.addr i
v2.addr j
x.addr i
w.addr i
x.addr i
A.addr ((j * 4096) + i)
x.addr i
z.addr i
A.addr ((i * 4096) + j)
x.addr j
w.addr i
BC Array cost info: Total number of arrays: 9
A.addr 9
u1.addr 2
u2.addr 2
v1.addr 2
v2.addr 2
w.addr 5
x.addr 12
y.addr 2
z.addr 2
BC Average cost: 4.222222e+00

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %n
double %alpha
double %beta
double* %A
double* %u1
double* %v1
double* %u2
double* %v2
double* %w
double* %x
double* %y
double* %z

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 4096)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 4096)
----Loop inc: (j + 1)
----Loop predicate: <
------array access A.addr ((i * 4096) + j)
------array access u1.addr i
------array access v1.addr j
------array access u2.addr i
------array access v2.addr j
------array access A.addr ((i * 4096) + j)
--i
--Loop Bound: (0, 4096)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 4096)
----Loop inc: (j + 1)
----Loop predicate: <
------array access x.addr i
------array access A.addr ((j * 4096) + i)
------array access y.addr j
------array access x.addr i
--i
--Loop Bound: (0, 4096)
--Loop inc: (i + 1)
--Loop predicate: <
----array access x.addr i
----array access z.addr i
----array access x.addr i
--i
--Loop Bound: (0, 4096)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 4096)
----Loop inc: (j + 1)
----Loop predicate: <
------array access w.addr i
------array access A.addr ((i * 4096) + j)
------array access x.addr j
------array access w.addr i

Finish analysis loops */ 
/* # of Out-most Loops: 4 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 40
------Sample number: 1677
----Sample number: 40
------Sample number: 1677
----Sample number: 40
----Sample number: 40
------Sample number: 1677
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* A_addr5	24576 */
/* A_addr0	24576 */
/* u1_addr1	24576 */
/* y_addr8	16384 */
/* x_addr9	16384 */
/* v1_addr2	24576 */
/* u2_addr3	24576 */
/* v2_addr4	24576 */
/* x_addr12	3 */
/* w_addr13	16384 */
/* x_addr6	16384 */
/* A_addr7	16384 */
/* x_addr10	3 */
/* z_addr11	3 */
/* A_addr14	16384 */
/* x_addr15	16384 */
/* w_addr16	16384 */
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
/* Array A_addr	i j */ 
/* i */
/* A_addr ((i * 4096) + j) 0 */
int calAddrA_addr0( int i, int j) {
    int result = (((i * 4096) + j)) * 8 / 64;
    return result;
}
/* Array u1_addr	i */ 
/* i */
/* u1_addr i 1 */
int calAddru1_addr1( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array v1_addr	j */ 
/* i */
/* v1_addr j 2 */
int calAddrv1_addr2( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* Array u2_addr	i */ 
/* i */
/* u2_addr i 3 */
int calAddru2_addr3( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array v2_addr	j */ 
/* i */
/* v2_addr j 4 */
int calAddrv2_addr4( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* Array A_addr	i j */ 
/* i */
/* A_addr ((i * 4096) + j) 5 */
int calAddrA_addr5( int i, int j) {
    int result = (((i * 4096) + j)) * 8 / 64;
    return result;
}
/* Array x_addr	i */ 
/* i */
/* x_addr i 6 */
int calAddrx_addr6( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array A_addr	j i */ 
/* i */
/* i */
/* i, j */
/* A_addr ((j * 4096) + i) 7 */
int calAddrA_addr7( int i, int j) {
    int result = (((j * 4096) + i)) * 8 / 64;
    return result;
}
/* Array y_addr	j */ 
/* i */
/* y_addr j 8 */
int calAddry_addr8( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* Array x_addr	i */ 
/* i */
/* x_addr i 9 */
int calAddrx_addr9( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array x_addr	i */ 
/* i */
/* x_addr i 10 */
int calAddrx_addr10( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array z_addr	i */ 
/* i */
/* z_addr i 11 */
int calAddrz_addr11( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array x_addr	i */ 
/* i */
/* x_addr i 12 */
int calAddrx_addr12( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array w_addr	i */ 
/* i */
/* w_addr i 13 */
int calAddrw_addr13( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array A_addr	i j */ 
/* i */
/* A_addr ((i * 4096) + j) 14 */
int calAddrA_addr14( int i, int j) {
    int result = (((i * 4096) + j)) * 8 / 64;
    return result;
}
/* Array x_addr	j */ 
/* i */
/* x_addr j 15 */
int calAddrx_addr15( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* Array w_addr	i */ 
/* i */
/* w_addr i 16 */
int calAddrw_addr16( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
void ref_A_addr5(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 4096; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j);
                            /* i 24576 */
                            /* j 6 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr5 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j);
                            /* i 24576 */
                            /* j 6 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr5 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        {
        int iLB2 = 0;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(4096, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            } else {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 16384, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr5 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr14( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr14, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(4096, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            } else {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            }
                            middle_accesses += 67121152;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 16384, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr5 --> A_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr14;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
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
void ref_A_addr0(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 4096; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr0, _1, j);
                            /* i 24576 */
                            /* j 6 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr0 --> A_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr0;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
                    if ( calAddrA_addr5( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr5, _1, j);
                            /* i 24576 */
                            /* j 6 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr0 --> A_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB2 = 0;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(4096, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            } else {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 16384, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr0 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr14( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr0, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr14, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(4096, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            } else {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            }
                            middle_accesses += 67121152;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 16384, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr0 --> A_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr14;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
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
void ref_u1_addr1(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 4096; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru1_addr1( i, j) == calAddru1_addr1(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru1_addr1, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru1_addr1, _1, j);
                            /* i 24576 */
                            /* j 6 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u1_addr1 --> u1_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u1_addr1;" << "(" << i_Start<< ", " << j_Start<< ");" << "u1_addr1;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB2 = 0;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
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
void ref_y_addr8(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr8( i, j) == calAddry_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry_addr8, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y_addr8 --> y_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "y_addr8;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
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
void ref_x_addr9(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr6( i, j) == calAddrx_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr6, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[x_addr9 --> x_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "x_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr6;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr9( i, j) == calAddrx_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr9, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[x_addr9 --> x_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "x_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr9;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr10( i) == calAddrx_addr9(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr9, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr10, _1);
                        /* i 3 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(4096, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 3;
                        } else {
                                    middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 3;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 3, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[x_addr9 --> x_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "x_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr10;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr12( i) == calAddrx_addr9(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr9, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr12, _1);
                        /* i 3 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(4096, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 3;
                        } else {
                                    middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 3;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 3, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[x_addr9 --> x_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "x_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr12;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr15( i, j) == calAddrx_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr15, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(4096, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            } else {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            }
                            middle_accesses += 12288;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[x_addr9 --> x_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "x_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr15;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_v1_addr2(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 4096; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv1_addr2( i, j) == calAddrv1_addr2(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv1_addr2, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv1_addr2, _1, j);
                            /* i 24576 */
                            /* j 6 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v1_addr2 --> v1_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v1_addr2;" << "(" << i_Start<< ", " << j_Start<< ");" << "v1_addr2;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int iLB2 = 0;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
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
void ref_u2_addr3(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 4096; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru2_addr3( i, j) == calAddru2_addr3(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru2_addr3, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru2_addr3, _1, j);
                            /* i 24576 */
                            /* j 6 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[u2_addr3 --> u2_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "u2_addr3;" << "(" << i_Start<< ", " << j_Start<< ");" << "u2_addr3;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        {
        int iLB2 = 0;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
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
void ref_v2_addr4(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 4096; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv2_addr4( i, j) == calAddrv2_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv2_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv2_addr4, _1, j);
                            /* i 24576 */
                            /* j 6 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[v2_addr4 --> v2_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "v2_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "v2_addr4;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int iLB2 = 0;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
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
void ref_x_addr12(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 40;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr10( i) == calAddrx_addr12(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr12, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr10, _1);
                        /* i 3 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 3, 3, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[x_addr12 --> x_addr10] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "x_addr12;" << "(" << i_Start<< ");" << "x_addr10;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr12( i) == calAddrx_addr12(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr12, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr12, _1);
                        /* i 3 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 3, 3, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[x_addr12 --> x_addr12] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "x_addr12;" << "(" << i_Start<< ");" << "x_addr12;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr15( i, j) == calAddrx_addr12(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr15, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(4096, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 3 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            } else {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 3 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 3, 16384, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[x_addr12 --> x_addr15] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "x_addr12;" << "(" << i_Start<< ");" << "x_addr15;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_w_addr13(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr13( i, j) == calAddrw_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrw_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrw_addr13, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[w_addr13 --> w_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "w_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "w_addr13;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr16( i, j) == calAddrw_addr13(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrw_addr13, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrw_addr16, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[w_addr13 --> w_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "w_addr13;" << "(" << i_Start<< ", " << j_Start<< ");" << "w_addr16;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_x_addr6(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr6( i, j) == calAddrx_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr6, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[x_addr6 --> x_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "x_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr6;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr9( i, j) == calAddrx_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr9, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[x_addr6 --> x_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "x_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr9;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr10( i) == calAddrx_addr6(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr6, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr10, _1);
                        /* i 3 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(4096, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 3;
                        } else {
                                    middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 3;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 3, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[x_addr6 --> x_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "x_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr10;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr12( i) == calAddrx_addr6(i_Start, j_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr6, _1, j_Start);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr12, _1);
                        /* i 3 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(4096, thread_cnt)) {
                                    middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 3;
                        } else {
                                    middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 3;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 3, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[x_addr6 --> x_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "x_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr12;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr15( i, j) == calAddrx_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr15, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(4096, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            } else {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            }
                            middle_accesses += 12288;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[x_addr6 --> x_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "x_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr15;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_A_addr7(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 4096; i=(i + 1)) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, j) == calAddrA_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr7, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr7 --> A_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr7;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr14( i, j) == calAddrA_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr14, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(4096, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            } else {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 16384 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            }
                            middle_accesses += 12288;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr7 --> A_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr14;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
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
void ref_x_addr10(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 40;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr10( i) == calAddrx_addr10(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr10, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr10, _1);
                        /* i 3 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 3, 3, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[x_addr10 --> x_addr10] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "x_addr10;" << "(" << i_Start<< ");" << "x_addr10;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr12( i) == calAddrx_addr10(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr10, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr12, _1);
                        /* i 3 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 3, 3, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[x_addr10 --> x_addr12] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "x_addr10;" << "(" << i_Start<< ");" << "x_addr12;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr15( i, j) == calAddrx_addr10(i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr10, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr15, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(4096, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 3 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            } else {
                                            middle_accesses += ((getChunkNum(4096, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (4096 % (CHUNK_SIZE * thread_cnt))) * 3 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 16384;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 3, 16384, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[x_addr10 --> x_addr15] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "x_addr10;" << "(" << i_Start<< ");" << "x_addr15;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_z_addr11(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 40;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        string idx_string =  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 4096; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrz_addr11( i) == calAddrz_addr11(i_Start)) {
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrz_addr11, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrz_addr11, _1);
                        /* i 3 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 3, 3, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                        rtHistoCal(RT, parallel_rt, 1.0);
#else
                        rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                            cout << "[z_addr11 --> z_addr11] [" << parallel_rt << "] (" << i_Start<< ") --> (" << i<< ") " << endl;
#endif
                            // cout << "z_addr11;" << "(" << i_Start<< ");" << "z_addr11;" << "(" << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
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
void ref_A_addr14(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr14( i, j) == calAddrA_addr14(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr14, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr14, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[A_addr14 --> A_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "A_addr14;" << "(" << i_Start<< ", " << j_Start<< ");" << "A_addr14;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
EndSample:
        s++;
        }
}
void ref_x_addr15(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr15( i, j) == calAddrx_addr15(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrx_addr15, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrx_addr15, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[x_addr15 --> x_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "x_addr15;" << "(" << i_Start<< ", " << j_Start<< ");" << "x_addr15;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_w_addr16(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1677;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 4096) { goto SAMPLE; }
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 4096; i=(i + 1)) {
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 4096; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr13( i, j) == calAddrw_addr16(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrw_addr16, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrw_addr13, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[w_addr16 --> w_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "w_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "w_addr13;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr16( i, j) == calAddrw_addr16(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrw_addr16, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrw_addr16, _1, j);
                            /* i 16384 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 16384, 16384, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[w_addr16 --> w_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "w_addr16;" << "(" << i_Start<< ", " << j_Start<< ");" << "w_addr16;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
    ref_A_addr5(thread_cnt, tmpRT);
    ref_A_addr0(thread_cnt, tmpRT);
    ref_u1_addr1(thread_cnt, tmpRT);
    ref_y_addr8(thread_cnt, tmpRT);
    ref_x_addr9(thread_cnt, tmpRT);
    ref_v1_addr2(thread_cnt, tmpRT);
    ref_u2_addr3(thread_cnt, tmpRT);
    ref_v2_addr4(thread_cnt, tmpRT);
    ref_x_addr12(thread_cnt, tmpRT);
    ref_w_addr13(thread_cnt, tmpRT);
    ref_x_addr6(thread_cnt, tmpRT);
    ref_A_addr7(thread_cnt, tmpRT);
    ref_x_addr10(thread_cnt, tmpRT);
    ref_z_addr11(thread_cnt, tmpRT);
    ref_A_addr14(thread_cnt, tmpRT);
    ref_x_addr15(thread_cnt, tmpRT);
    ref_w_addr16(thread_cnt, tmpRT);
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
 /* Analyze function: gemver */ 
