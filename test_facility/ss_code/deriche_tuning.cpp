
 /* Start to analysis array index
Array index info: Total number of references: 20
y1.addr ((i * 8192) + j)
imgIn.addr ((i * 8192) + j)
y1.addr ((i * 8192) + j)
y2.addr ((((i * 8192) + 8192) - 1) - j)
imgIn.addr ((i * 8192) + j)
y1.addr ((i * 8192) + j)
y2.addr ((i * 8192) + j)
y1.addr ((i * 8192) + j)
imgOut.addr ((i * 8192) + j)
y1.addr ((i * 8192) + j)
y2.addr ((((i * 8192) + 8192) - 1) - j)
imgIn.addr ((((i * 8192) + 8192) - 1) - j)
imgOut.addr ((i * 8192) + j)
imgOut.addr ((i * 8192) + j)
y2.addr (((8191 - i) * 8192) + j)
imgOut.addr (((8191 - i) * 8192) + j)
y2.addr (((8191 - i) * 8192) + j)
y1.addr ((i * 8192) + j)
y2.addr ((i * 8192) + j)
imgOut.addr ((i * 8192) + j)
BC Array cost info: Total number of arrays: 4
imgIn.addr 6
imgOut.addr 12
y1.addr 14
y2.addr 14
BC Average cost: 1.150000e+01

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %y1
double* %imgIn
double* %y2
double* %imgOut
double %alpha

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 8192)
----Loop inc: (j + 1)
----Loop predicate: <
------array access imgIn.addr ((i * 8192) + j)
------array access y1.addr ((i * 8192) + j)
------array access imgIn.addr ((i * 8192) + j)
------array access y1.addr ((i * 8192) + j)
--i
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 8192)
----Loop inc: (j + 1)
----Loop predicate: <
------array access y2.addr ((((i * 8192) + 8192) - 1) - j)
------array access imgIn.addr ((((i * 8192) + 8192) - 1) - j)
------array access y2.addr ((((i * 8192) + 8192) - 1) - j)
--i
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 8192)
----Loop inc: (j + 1)
----Loop predicate: <
------array access y1.addr ((i * 8192) + j)
------array access y2.addr ((i * 8192) + j)
------array access imgOut.addr ((i * 8192) + j)
--j
--Loop Bound: (0, 8192)
--Loop inc: (j + 1)
--Loop predicate: <
----i
----Loop Bound: (0, 8192)
----Loop inc: (i + 1)
----Loop predicate: <
------array access imgOut.addr ((i * 8192) + j)
------array access y1.addr ((i * 8192) + j)
------array access imgOut.addr ((i * 8192) + j)
------array access y1.addr ((i * 8192) + j)
--j
--Loop Bound: (0, 8192)
--Loop inc: (j + 1)
--Loop predicate: <
----i
----Loop Bound: (0, 8192)
----Loop inc: (i + 1)
----Loop predicate: <
------array access y2.addr (((8191 - i) * 8192) + j)
------array access imgOut.addr (((8191 - i) * 8192) + j)
------array access y2.addr (((8191 - i) * 8192) + j)
--i
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 8192)
----Loop inc: (j + 1)
----Loop predicate: <
------array access y1.addr ((i * 8192) + j)
------array access y2.addr ((i * 8192) + j)
------array access imgOut.addr ((i * 8192) + j)

Finish analysis loops */ 
/* # of Out-most Loops: 6 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 81
------Sample number: 6710
----Sample number: 81
------Sample number: 6710
----Sample number: 81
------Sample number: 6710
----Sample number: 81
------Sample number: 6710
----Sample number: 81
------Sample number: 6710
----Sample number: 81
------Sample number: 6710
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* y1_addr1	32768 */
/* imgIn_addr2	32768 */
/* y1_addr3	32768 */
/* y2_addr6	24576 */
/* imgIn_addr0	32768 */
/* y1_addr7	24576 */
/* y2_addr8	24576 */
/* y1_addr11	32768 */
/* imgOut_addr12	32768 */
/* y1_addr13	32768 */
/* y2_addr4	24576 */
/* imgIn_addr5	24576 */
/* imgOut_addr9	24576 */
/* imgOut_addr10	32768 */
/* y2_addr14	24576 */
/* imgOut_addr15	24576 */
/* y2_addr16	24576 */
/* y1_addr17	24576 */
/* y2_addr18	24576 */
/* imgOut_addr19	24576 */
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
/* Array imgIn_addr	i j */ 
/* i */
/* imgIn_addr ((i * 8192) + j) 0 */
int calAddrimgIn_addr0( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y1_addr	i j */ 
/* i */
/* y1_addr ((i * 8192) + j) 1 */
int calAddry1_addr1( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array imgIn_addr	i j */ 
/* i */
/* imgIn_addr ((i * 8192) + j) 2 */
int calAddrimgIn_addr2( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y1_addr	i j */ 
/* i */
/* y1_addr ((i * 8192) + j) 3 */
int calAddry1_addr3( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y2_addr	i j */ 
/* i */
/* y2_addr ((((i * 8192) + 8192) - 1) - j) 4 */
int calAddry2_addr4( int i, int j) {
    int result = (((((i * 8192) + 8192) - 1) - j)) * 8 / 64;
    return result;
}
/* Array imgIn_addr	i j */ 
/* i */
/* imgIn_addr ((((i * 8192) + 8192) - 1) - j) 5 */
int calAddrimgIn_addr5( int i, int j) {
    int result = (((((i * 8192) + 8192) - 1) - j)) * 8 / 64;
    return result;
}
/* Array y2_addr	i j */ 
/* i */
/* y2_addr ((((i * 8192) + 8192) - 1) - j) 6 */
int calAddry2_addr6( int i, int j) {
    int result = (((((i * 8192) + 8192) - 1) - j)) * 8 / 64;
    return result;
}
/* Array y1_addr	i j */ 
/* i */
/* y1_addr ((i * 8192) + j) 7 */
int calAddry1_addr7( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y2_addr	i j */ 
/* i */
/* y2_addr ((i * 8192) + j) 8 */
int calAddry2_addr8( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array imgOut_addr	i j */ 
/* i */
/* imgOut_addr ((i * 8192) + j) 9 */
int calAddrimgOut_addr9( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array imgOut_addr	i j */ 
/* j */
/* j */
/* j, i */
/* imgOut_addr ((i * 8192) + j) 10 */
int calAddrimgOut_addr10( int j, int i) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y1_addr	i j */ 
/* j */
/* j */
/* j, i */
/* y1_addr ((i * 8192) + j) 11 */
int calAddry1_addr11( int j, int i) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array imgOut_addr	i j */ 
/* j */
/* j */
/* j, i */
/* imgOut_addr ((i * 8192) + j) 12 */
int calAddrimgOut_addr12( int j, int i) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y1_addr	i j */ 
/* j */
/* j */
/* j, i */
/* y1_addr ((i * 8192) + j) 13 */
int calAddry1_addr13( int j, int i) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y2_addr	i j */ 
/* j */
/* j */
/* j, i */
/* y2_addr (((8191 - i) * 8192) + j) 14 */
int calAddry2_addr14( int j, int i) {
    int result = ((((8191 - i) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array imgOut_addr	i j */ 
/* j */
/* j */
/* j, i */
/* imgOut_addr (((8191 - i) * 8192) + j) 15 */
int calAddrimgOut_addr15( int j, int i) {
    int result = ((((8191 - i) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y2_addr	i j */ 
/* j */
/* j */
/* j, i */
/* y2_addr (((8191 - i) * 8192) + j) 16 */
int calAddry2_addr16( int j, int i) {
    int result = ((((8191 - i) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y1_addr	i j */ 
/* i */
/* y1_addr ((i * 8192) + j) 17 */
int calAddry1_addr17( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array y2_addr	i j */ 
/* i */
/* y2_addr ((i * 8192) + j) 18 */
int calAddry2_addr18( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array imgOut_addr	i j */ 
/* i */
/* imgOut_addr ((i * 8192) + j) 19 */
int calAddrimgOut_addr19( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
void ref_y1_addr1(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8192; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr1( i, j) == calAddry1_addr1(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr1, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr1, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr1 --> y1_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr1;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr1;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr3( i, j) == calAddry1_addr1(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr1, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr3, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr1 --> y1_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr1;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr3;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB2 = 0;
        for ( int i = iLB2; i < 8192; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr7( i, j) == calAddry1_addr1(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr1, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr7, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 201326592;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr1 --> y1_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr1;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr7;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr1(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr1, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr11, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            }
                            middle_accesses += 402653184;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr1 --> y1_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr1;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr11;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr1(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr1, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr13, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            }
                            middle_accesses += 402653184;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr1 --> y1_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr1;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr13;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr1(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr1, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr17, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 872415232;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr1 --> y1_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr1;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr17;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_imgIn_addr2(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8192; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr0( i, j) == calAddrimgIn_addr2(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgIn_addr2, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgIn_addr0, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgIn_addr2 --> imgIn_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgIn_addr2;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgIn_addr0;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr2( i, j) == calAddrimgIn_addr2(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgIn_addr2, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgIn_addr2, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgIn_addr2 --> imgIn_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgIn_addr2;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgIn_addr2;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        for ( int i = iLB2; i < 8192; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr5( i, j) == calAddrimgIn_addr2(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgIn_addr2, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgIn_addr5, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgIn_addr2 --> imgIn_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgIn_addr2;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgIn_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
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
void ref_y1_addr3(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8192; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr1( i, j) == calAddry1_addr3(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr3, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr1, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr3 --> y1_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr3;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr1;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr3( i, j) == calAddry1_addr3(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr3, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr3, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr3 --> y1_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr3;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr3;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        for ( int i = iLB2; i < 8192; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr7( i, j) == calAddry1_addr3(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr3, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr7, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 201326592;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr3 --> y1_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr3;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr7;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr3(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr3, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr11, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            }
                            middle_accesses += 402653184;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr3 --> y1_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr3;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr11;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr3(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr3, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr13, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            }
                            middle_accesses += 402653184;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr3 --> y1_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr3;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr13;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr3(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr3, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr17, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 872415232;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr3 --> y1_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr3;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr17;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_y2_addr6(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 8192; i=(i + 1)) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr4( i, j) == calAddry2_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr4, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[y2_addr6 --> y2_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr4;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr6( i, j) == calAddry2_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr6, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[y2_addr6 --> y2_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr6;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr8( i, j) == calAddry2_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr8, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr6 --> y2_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr8;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr14, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 469762048;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr6 --> y2_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr14;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr16, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 469762048;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr6 --> y2_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr16;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr6(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr18, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 671088640;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr6 --> y2_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr6;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr18;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_imgIn_addr0(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8192; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr0( i, j) == calAddrimgIn_addr0(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgIn_addr0, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgIn_addr0, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgIn_addr0 --> imgIn_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgIn_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgIn_addr0;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr2( i, j) == calAddrimgIn_addr0(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgIn_addr0, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgIn_addr2, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgIn_addr0 --> imgIn_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgIn_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgIn_addr2;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB2 = 0;
        for ( int i = iLB2; i < 8192; i=(i + 1)) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr5( i, j) == calAddrimgIn_addr0(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgIn_addr0, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgIn_addr5, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgIn_addr0 --> imgIn_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgIn_addr0;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgIn_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
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
void ref_y1_addr7(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr7( i, j) == calAddry1_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr7, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[y1_addr7 --> y1_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr7;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr11, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 32768, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr7 --> y1_addr11] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr11;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr13, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 32768, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr7 --> y1_addr13] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr13;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr7(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr7, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr17, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 469762048;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr7 --> y1_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr7;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr17;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_y2_addr8(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr8( i, j) == calAddry2_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr8, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[y2_addr8 --> y2_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr8;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr14, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 268435456;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr8 --> y2_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr14;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr16, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 268435456;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr8 --> y2_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr16;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr8(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr8, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr18, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 469762048;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr8 --> y2_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr8;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr18;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_y1_addr11(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if (j_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(j_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB6 = j_Start;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr11(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr11, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr11, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr11 --> y1_addr11] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr11;" << "(" << j_Start<< ", " << i_Start<< ");" << "y1_addr11;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr11(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr11, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr13, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr11 --> y1_addr13] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr11;" << "(" << j_Start<< ", " << i_Start<< ");" << "y1_addr13;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr11(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr11, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr17, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 201326592;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr11 --> y1_addr17] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr11;" << "(" << j_Start<< ", " << i_Start<< ");" << "y1_addr17;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_imgOut_addr12(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if (j_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(j_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB6 = j_Start;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr10( j, i) == calAddrimgOut_addr12(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr12, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr10, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr12 --> imgOut_addr10] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr12;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr10;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr12( j, i) == calAddrimgOut_addr12(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr12, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr12, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr12 --> imgOut_addr12] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr12;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr12;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr15( j, i) == calAddrimgOut_addr12(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr12, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr15, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr12 --> imgOut_addr15] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr12;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr15;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr12(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr12, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr19, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 201326592;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr12 --> imgOut_addr19] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgOut_addr12;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr19;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_y1_addr13(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if (j_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(j_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB6 = j_Start;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr13(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr13, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr11, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr13 --> y1_addr11] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr13;" << "(" << j_Start<< ", " << i_Start<< ");" << "y1_addr11;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr13(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr13, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr13, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr13 --> y1_addr13] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y1_addr13;" << "(" << j_Start<< ", " << i_Start<< ");" << "y1_addr13;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr13(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr13, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr17, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 201326592;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y1_addr13 --> y1_addr17] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr13;" << "(" << j_Start<< ", " << i_Start<< ");" << "y1_addr17;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_y2_addr4(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 8192; i=(i + 1)) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr4( i, j) == calAddry2_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr4, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[y2_addr4 --> y2_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr4;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr6( i, j) == calAddry2_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr6, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[y2_addr4 --> y2_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr6;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr8( i, j) == calAddry2_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr8, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr4 --> y2_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr8;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr14, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 469762048;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr4 --> y2_addr14] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr14;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr16, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 469762048;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr4 --> y2_addr16] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr16;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr4(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr4, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr18, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 671088640;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr4 --> y2_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr4;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr18;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_imgIn_addr5(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 8192; i=(i + 1)) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr5( i, j) == calAddrimgIn_addr5(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgIn_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgIn_addr5, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[imgIn_addr5 --> imgIn_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgIn_addr5;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgIn_addr5;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
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
void ref_imgOut_addr9(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 8192; i=(i + 1)) {
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr9( i, j) == calAddrimgOut_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr9, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[imgOut_addr9 --> imgOut_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgOut_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgOut_addr9;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        {
        int jLB6 = 0;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr10( j, i) == calAddrimgOut_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr10, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 32768, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr9 --> imgOut_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgOut_addr10;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr12( j, i) == calAddrimgOut_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr12, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 32768;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 32768, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr9 --> imgOut_addr12] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgOut_addr12;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr15( j, i) == calAddrimgOut_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr15, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 268435456;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr9 --> imgOut_addr15] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgOut_addr15;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr9(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr9, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr19, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((i_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 469762048;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr9 --> imgOut_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgOut_addr9;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgOut_addr19;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_imgOut_addr10(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if (j_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(j_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB6 = j_Start;
        for ( int j = jLB6; j < 8192; j=(j + 1)) {
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 8192; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr10( j, i) == calAddrimgOut_addr10(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr10, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr10, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr10 --> imgOut_addr10] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr10;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr10;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr12( j, i) == calAddrimgOut_addr10(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr10, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr12, _1, i);
                            /* j 32768 */
                            /* i 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 32768, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr10 --> imgOut_addr12] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr10;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr12;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr15( j, i) == calAddrimgOut_addr10(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr10, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr15, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((j - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr10 --> imgOut_addr15] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr10;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr15;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr10(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr10, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr19, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 32768 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 201326592;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (i - 0), cnt, 32768, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr10 --> imgOut_addr19] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgOut_addr10;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr19;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_y2_addr14(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if (j_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(j_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB8 = j_Start;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            if ( j == j_Start ) {
                iLB9 = i_Start;
            }
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr14(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr14, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr14, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr14 --> y2_addr14] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr14;" << "(" << j_Start<< ", " << i_Start<< ");" << "y2_addr14;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr14(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr14, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr16, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr14 --> y2_addr16] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr14;" << "(" << j_Start<< ", " << i_Start<< ");" << "y2_addr16;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr14(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr14, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr18, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr14 --> y2_addr18] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr14;" << "(" << j_Start<< ", " << i_Start<< ");" << "y2_addr18;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_imgOut_addr15(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if (j_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(j_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB8 = j_Start;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            if ( j == j_Start ) {
                iLB9 = i_Start;
            }
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr15( j, i) == calAddrimgOut_addr15(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr15, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr15, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr15 --> imgOut_addr15] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "imgOut_addr15;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr15;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr15(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr15, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr19, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[imgOut_addr15 --> imgOut_addr19] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgOut_addr15;" << "(" << j_Start<< ", " << i_Start<< ");" << "imgOut_addr19;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_y2_addr16(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        if (j_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(j_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB8 = j_Start;
        for ( int j = jLB8; j < 8192; j=(j + 1)) {
            {
            int iLB9 = 0;
            if ( j == j_Start ) {
                iLB9 = i_Start;
            }
            for ( int i = iLB9; i < 8192; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr16(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr16, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr14, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr16 --> y2_addr14] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr16;" << "(" << j_Start<< ", " << i_Start<< ");" << "y2_addr14;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr16(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr16, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr16, _1, i);
                            /* j 24576 */
                            /* i 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (j - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr16 --> y2_addr16] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << j<< ", " << i<< ") " << endl;
#endif
                                // cout << "y2_addr16;" << "(" << j_Start<< ", " << i_Start<< ");" << "y2_addr16;" << "(" << j<< ", " << i<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr16(j_Start, i_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr16, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr18, _1, j);
                            /* i 24576 */
                            /* j 3 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192, thread_cnt)) {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            } else {
                                            middle_accesses += ((getChunkNum(8192, thread_cnt) - getChunkID((j_Start -0), thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (8192 % (CHUNK_SIZE * thread_cnt))) * 24576 + getChunkID((i - 0), thread_cnt) * CHUNK_SIZE * thread_cnt * 24576;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            uint64_t parallel_rt = parallel_predict((j_Start -0), (i - 0), cnt, 24576, 24576, middle_accesses, thread_cnt, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
#ifdef DEBUG
                            rtHistoCal(RT, parallel_rt, 1.0);
#else
                            rtHistoCal(RT, parallel_rt, 1.0);
#endif
#ifdef DEBUG
                                cout << "[y2_addr16 --> y2_addr18] [" << parallel_rt << "] (" << j_Start<< ", " << i_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr16;" << "(" << j_Start<< ", " << i_Start<< ");" << "y2_addr18;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_y1_addr17(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB10 = i_Start;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            if ( i == i_Start ) {
                jLB11 = j_Start;
            }
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr17(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry1_addr17, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry1_addr17, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[y1_addr17 --> y1_addr17] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y1_addr17;" << "(" << i_Start<< ", " << j_Start<< ");" << "y1_addr17;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_y2_addr18(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB10 = i_Start;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            if ( i == i_Start ) {
                jLB11 = j_Start;
            }
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr18(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddry2_addr18, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddry2_addr18, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[y2_addr18 --> y2_addr18] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "y2_addr18;" << "(" << i_Start<< ", " << j_Start<< ");" << "y2_addr18;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
void ref_imgOut_addr19(int thread_cnt, map<uint64_t, double> & RT) {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6710;) {
SAMPLE:
        int i_Start = rand() % (8192 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + thread_cnt * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (8192 - 0) + 0;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB10 = i_Start;
        for ( int i = iLB10; i < 8192; i=(i + 1)) {
            {
            int jLB11 = 0;
            if ( i == i_Start ) {
                jLB11 = j_Start;
            }
            for ( int j = jLB11; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr19(i_Start, j_Start)) {
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrimgOut_addr19, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrimgOut_addr19, _1, j);
                            /* i 24576 */
                            /* j 3 */
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
                                cout << "[imgOut_addr19 --> imgOut_addr19] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                                // cout << "imgOut_addr19;" << "(" << i_Start<< ", " << j_Start<< ");" << "imgOut_addr19;" << "(" << i<< ", " << j<< ");" << cnt << ";" << parallel_rt << ";" << reuse_type << endl;
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
    ref_y1_addr1(thread_cnt, tmpRT);
    ref_imgIn_addr2(thread_cnt, tmpRT);
    ref_y1_addr3(thread_cnt, tmpRT);
    ref_y2_addr6(thread_cnt, tmpRT);
    ref_imgIn_addr0(thread_cnt, tmpRT);
    ref_y1_addr7(thread_cnt, tmpRT);
    ref_y2_addr8(thread_cnt, tmpRT);
    ref_y1_addr11(thread_cnt, tmpRT);
    ref_imgOut_addr12(thread_cnt, tmpRT);
    ref_y1_addr13(thread_cnt, tmpRT);
    ref_y2_addr4(thread_cnt, tmpRT);
    ref_imgIn_addr5(thread_cnt, tmpRT);
    ref_imgOut_addr9(thread_cnt, tmpRT);
    ref_imgOut_addr10(thread_cnt, tmpRT);
    ref_y2_addr14(thread_cnt, tmpRT);
    ref_imgOut_addr15(thread_cnt, tmpRT);
    ref_y2_addr16(thread_cnt, tmpRT);
    ref_y1_addr17(thread_cnt, tmpRT);
    ref_y2_addr18(thread_cnt, tmpRT);
    ref_imgOut_addr19(thread_cnt, tmpRT);
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
 /* Analyze function: deriche */ 