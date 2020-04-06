
 /* Start to analysis array index
Array index info: Total number of references: 34
p.addr (((i * 1024) + j) - 1)
u.addr ((j * 1024) + i)
u.addr (((j * 1024) + i) + 1)
v.addr (0 + i)
p.addr ((i * 1024) + 0)
v.addr (0 + i)
q.addr ((i * 1024) + 0)
v.addr (1047552 + i)
p.addr ((i * 1024) + j)
v.addr (((j + 1) * 1024) + i)
q.addr ((i * 1024) + j)
p.addr ((i * 1024) + j)
u.addr ((i * 1024) + 0)
u.addr (((j * 1024) + i) - 1)
q.addr ((i * 1024) + 0)
p.addr (((i * 1024) + j) - 1)
p.addr ((i * 1024) + j)
v.addr (((i - 1) * 1024) + j)
q.addr (((i * 1024) + j) - 1)
p.addr (((i * 1024) + j) - 1)
q.addr ((i * 1024) + j)
u.addr (((i * 1024) + j) + 1)
q.addr ((i * 1024) + j)
u.addr ((i * 1024) + j)
v.addr ((j * 1024) + i)
p.addr ((i * 1024) + 0)
u.addr ((i * 1024) + 0)
v.addr ((i * 1024) + j)
v.addr (((i + 1) * 1024) + j)
q.addr (((i * 1024) + j) - 1)
p.addr (((i * 1024) + j) - 1)
q.addr ((i * 1024) + j)
u.addr (((i * 1024) + 1024) - 1)
p.addr ((i * 1024) + j)

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
--t
--Loop Bound: (1, 10)
--Loop inc: (t + 1)
--Loop predicate: <=
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------array access v.addr (0 + i)
------array access p.addr ((i * 1024) + 0)
------array access v.addr (0 + i)
------array access q.addr ((i * 1024) + 0)
------j
------Loop Bound: (1, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access p.addr (((i * 1024) + j) - 1)
--------array access p.addr ((i * 1024) + j)
--------array access u.addr (((j * 1024) + i) - 1)
--------array access u.addr ((j * 1024) + i)
--------array access u.addr (((j * 1024) + i) + 1)
--------array access q.addr (((i * 1024) + j) - 1)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access q.addr ((i * 1024) + j)
------array access v.addr (1047552 + i)
------j
------Loop Bound: (1022, 1)
------Loop inc: (j + -1)
------Loop predicate: >=
--------array access p.addr ((i * 1024) + j)
--------array access v.addr (((j + 1) * 1024) + i)
--------array access q.addr ((i * 1024) + j)
--------array access v.addr ((j * 1024) + i)
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------array access u.addr ((i * 1024) + 0)
------array access p.addr ((i * 1024) + 0)
------array access u.addr ((i * 1024) + 0)
------array access q.addr ((i * 1024) + 0)
------j
------Loop Bound: (1, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access p.addr (((i * 1024) + j) - 1)
--------array access p.addr ((i * 1024) + j)
--------array access v.addr (((i - 1) * 1024) + j)
--------array access v.addr ((i * 1024) + j)
--------array access v.addr (((i + 1) * 1024) + j)
--------array access q.addr (((i * 1024) + j) - 1)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access q.addr ((i * 1024) + j)
------array access u.addr (((i * 1024) + 1024) - 1)
------j
------Loop Bound: (1022, 1)
------Loop inc: (j + -1)
------Loop predicate: >=
--------array access p.addr ((i * 1024) + j)
--------array access u.addr (((i * 1024) + j) + 1)
--------array access q.addr ((i * 1024) + j)
--------array access u.addr ((i * 1024) + j)

Finish analysis loops */ 
/* # of Out-most Loops: 1 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 1
------Sample number: 102
--------Sample number: 10444
--------Sample number: 10444
------Sample number: 102
--------Sample number: 10444
--------Sample number: 10444
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* p_addr4	16721964 */
/* u_addr7	16721964 */
/* u_addr8	16721964 */
/* v_addr0	16721964 */
/* p_addr1	16721964 */
/* v_addr2	16721964 */
/* q_addr3	16721964 */
/* v_addr12	16721964 */
/* p_addr13	16721964 */
/* v_addr14	16721964 */
/* q_addr15	16721964 */
/* p_addr5	16721964 */
/* u_addr17	16721964 */
/* u_addr6	16721964 */
/* q_addr20	16721964 */
/* p_addr21	16721964 */
/* p_addr22	16721964 */
/* v_addr23	16721964 */
/* q_addr9	16721964 */
/* p_addr10	16721964 */
/* q_addr11	16721964 */
/* u_addr31	16721964 */
/* q_addr32	16721964 */
/* u_addr33	16721964 */
/* v_addr16	16721964 */
/* p_addr18	16721964 */
/* u_addr19	16721964 */
/* v_addr24	16721964 */
/* v_addr25	16721964 */
/* q_addr26	16721964 */
/* p_addr27	16721964 */
/* q_addr28	16721964 */
/* u_addr29	16721964 */
/* p_addr30	16721964 */
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <functional>
#ifdef PAPI_TIMER
#  include <chrono>
#endif
using namespace std;
using namespace placeholders;
#ifdef PAPI_TIMER
using namespace  chrono;
#endif
 map<uint64_t, double> RT;
 map<uint64_t, double> MR;
int getChunkNum(uint64_t is) {
    return is / (CHUNK_SIZE * THREAD_NUM);
}
int getChunkID(uint64_t i) {
    return floor(i / (CHUNK_SIZE * THREAD_NUM));
}
int getThreadID(uint64_t i) {
    return i / CHUNK_SIZE - floor(i / (CHUNK_SIZE * THREAD_NUM))*THREAD_NUM ;
}
int getThreadLocalPos(uint64_t i) {
    return i % CHUNK_SIZE;
}
int search_src_candidate_neighbor(uint64_t i, function<uint64_t(uint64_t)> calAddr) {
    int c_start = i % CHUNK_SIZE + i / (CHUNK_SIZE * THREAD_NUM) * THREAD_NUM * CHUNK_SIZE;
    int c_end = c_start + (THREAD_NUM - 1) * CHUNK_SIZE;
    for (int c = i + CHUNK_SIZE; c <= c_end; c=c+CHUNK_SIZE) {
        if (calAddr(i) == calAddr(c)) { return getThreadID(c); }
    }
    return -1;
}
int search_sink_candidate_neighbor(uint64_t i, function<uint64_t(uint64_t)> calAddr) {
    int c_start = i % CHUNK_SIZE + i / (CHUNK_SIZE * THREAD_NUM) * THREAD_NUM * CHUNK_SIZE;
    for (int c = c_start; c <= i; c=c+CHUNK_SIZE) {
        if (calAddr(i) == calAddr(c)) { return getThreadID(c); }
    }
    return -1;
}
uint64_t parallel_predict(uint64_t i_src, uint64_t i_sink, uint64_t rt, uint64_t lsrc, uint64_t lsink, uint64_t middle_accesses, bool is_normal_ref, bool is_in_same_loop, function<uint64_t(uint64_t)> srcAddrCal, function<uint64_t(uint64_t)> sinkAddrCal) {
    uint64_t parallel_rt = rt;
    int tsrc = getThreadID(i_src);
    int tsink = getThreadID(i_sink);
    int dT = tsink - tsrc;
    if (!is_in_same_loop || getChunkID(i_src) != getChunkID(i_sink)) {
#ifdef DEBUG
        cout << "Inter Chunk Reuse" << endl;
#endif
        cout << "rt " << rt << endl;
#endif
        parallel_rt = rt * THREAD_NUM - CHUNK_SIZE * THREAD_NUM * (lsrc*(THREAD_NUM - tsrc) + lsink * tsink) + CHUNK_SIZE * THREAD_NUM * lsrc - (THREAD_NUM - 1) * middle_accesses + dT;
    } else if (!is_normal_ref) {
        /* intra chunk reuse */
#ifdef DEBUG
        cout << "Neighboring Effect" << endl;
#endif
        int tsrc_neighbor = search_src_candidate_neighbor(i_src, srcAddrCal);
        int tsink_neighbor = search_sink_candidate_neighbor(i_sink, sinkAddrCal);
        if (tsrc_neighbor >= 0) {
#ifdef DEBUG
            cout << "Find sink in src neighbor at" << tsrc_neighbor << endl;
#endif
            return tsrc_neighbor - tsrc;
        } else if (tsink_neighbor >= 0) {
#ifdef DEBUG
            cout << "Find sink in sink neighbor at" << tsink_neighbor << endl;
#endif
            if (getChunkID(i_src) == getChunkID(i_sink)) { return rt * THREAD_NUM + tsink_neighbor - tsink; }
        }
    } else if (getChunkID(i_src) == getChunkID(i_sink)) {
        /* same thread -- scaling effect */
        if (dT == 0) {
#ifdef DEBUG
            cout << "Scaling Effect" << endl;
#endif
            parallel_rt = rt * THREAD_NUM;
        } else if (getThreadLocalPos(i_src) <= getThreadLocalPos(i_sink)) { // src-sink order
            if ((rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + dT) < 0) { printf("NORMAL ORDER NEGATIVE PRI\n"); }
#ifdef DEBUG
            cout << "Src-Sink Order Folding Effect" << endl;
#endif
            parallel_rt = rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + abs(dT);
        } else { // sink-src order
            if ((rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + dT) < 0) { printf("REVERSE ORDER NEGATIVE PRI\n"); }
#ifdef DEBUG
            cout << "Sink-Src Order Folding Effect" << endl;
#endif
            // parallel_rt = CHUNK_SIZE * lsrc * THREAD_NUM * dT - (rt * THREAD_NUM) - abs(dT);
            return 0;
        }
    }
    return parallel_rt;
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
/* Array v_addr	i */ 
/* v_addr (0 + i) 0 */
int calAddrv_addr0( int t, int i) {
    int result = ((0 + i)) * 8 / 64;
    return result;
}
/* Array p_addr	i */ 
/* p_addr ((i * 1024) + 0) 1 */
int calAddrp_addr1( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* Array v_addr	i */ 
/* v_addr (0 + i) 2 */
int calAddrv_addr2( int t, int i) {
    int result = ((0 + i)) * 8 / 64;
    return result;
}
/* Array q_addr	i */ 
/* q_addr ((i * 1024) + 0) 3 */
int calAddrq_addr3( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* p_addr (((i * 1024) + j) - 1) 4 */
int calAddrp_addr4( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* p_addr ((i * 1024) + j) 5 */
int calAddrp_addr5( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array u_addr	j i */ 
/* u_addr (((j * 1024) + i) - 1) 6 */
int calAddru_addr6( int t, int i, int j) {
    int result = ((((j * 1024) + i) - 1)) * 8 / 64;
    return result;
}
/* Array u_addr	j i */ 
/* u_addr ((j * 1024) + i) 7 */
int calAddru_addr7( int t, int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 64;
    return result;
}
/* Array u_addr	j i */ 
/* u_addr (((j * 1024) + i) + 1) 8 */
int calAddru_addr8( int t, int i, int j) {
    int result = ((((j * 1024) + i) + 1)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* q_addr (((i * 1024) + j) - 1) 9 */
int calAddrq_addr9( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* p_addr (((i * 1024) + j) - 1) 10 */
int calAddrp_addr10( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* q_addr ((i * 1024) + j) 11 */
int calAddrq_addr11( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array v_addr	i */ 
/* v_addr (1047552 + i) 12 */
int calAddrv_addr12( int t, int i) {
    int result = ((1047552 + i)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* p_addr ((i * 1024) + j) 13 */
int calAddrp_addr13( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array v_addr	j i */ 
/* v_addr (((j + 1) * 1024) + i) 14 */
int calAddrv_addr14( int t, int i, int j) {
    int result = ((((j + 1) * 1024) + i)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* q_addr ((i * 1024) + j) 15 */
int calAddrq_addr15( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array v_addr	j i */ 
/* v_addr ((j * 1024) + i) 16 */
int calAddrv_addr16( int t, int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 64;
    return result;
}
/* Array u_addr	i */ 
/* u_addr ((i * 1024) + 0) 17 */
int calAddru_addr17( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* Array p_addr	i */ 
/* p_addr ((i * 1024) + 0) 18 */
int calAddrp_addr18( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* Array u_addr	i */ 
/* u_addr ((i * 1024) + 0) 19 */
int calAddru_addr19( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* Array q_addr	i */ 
/* q_addr ((i * 1024) + 0) 20 */
int calAddrq_addr20( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* p_addr (((i * 1024) + j) - 1) 21 */
int calAddrp_addr21( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* p_addr ((i * 1024) + j) 22 */
int calAddrp_addr22( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array v_addr	i j */ 
/* v_addr (((i - 1) * 1024) + j) 23 */
int calAddrv_addr23( int t, int i, int j) {
    int result = ((((i - 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* Array v_addr	i j */ 
/* v_addr ((i * 1024) + j) 24 */
int calAddrv_addr24( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array v_addr	i j */ 
/* v_addr (((i + 1) * 1024) + j) 25 */
int calAddrv_addr25( int t, int i, int j) {
    int result = ((((i + 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* q_addr (((i * 1024) + j) - 1) 26 */
int calAddrq_addr26( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* p_addr (((i * 1024) + j) - 1) 27 */
int calAddrp_addr27( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* q_addr ((i * 1024) + j) 28 */
int calAddrq_addr28( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array u_addr	i */ 
/* u_addr (((i * 1024) + 1024) - 1) 29 */
int calAddru_addr29( int t, int i) {
    int result = ((((i * 1024) + 1024) - 1)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* p_addr ((i * 1024) + j) 30 */
int calAddrp_addr30( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array u_addr	i j */ 
/* u_addr (((i * 1024) + j) + 1) 31 */
int calAddru_addr31( int t, int i, int j) {
    int result = ((((i * 1024) + j) + 1)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* q_addr ((i * 1024) + j) 32 */
int calAddrq_addr32( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array u_addr	i j */ 
/* u_addr ((i * 1024) + j) 33 */
int calAddru_addr33( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_p_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr4 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr4 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr4 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr4 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr4 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr4 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr4 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr4 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr4 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_u_addr7() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr7 --> u_addr6] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr7 --> u_addr7] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr7 --> u_addr8] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr7 --> u_addr17] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr7 --> u_addr19] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
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
                    if ( calAddru_addr29( t, i) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr7 --> u_addr29] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr7 --> u_addr31] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr7 --> u_addr33] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_u_addr8() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr8 --> u_addr6] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr8 --> u_addr7] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr8 --> u_addr8] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr8 --> u_addr17] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr8 --> u_addr19] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
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
                    if ( calAddru_addr29( t, i) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr8 --> u_addr29] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr8 --> u_addr31] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr8 --> u_addr33] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_v_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr0(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr0 --> v_addr0] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr0(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr0 --> v_addr2] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
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
                    if ( calAddrv_addr12( t, i) == calAddrv_addr0(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr0 --> v_addr12] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr0 --> v_addr14] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr0 --> v_addr16] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr0 --> v_addr23] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr0 --> v_addr24] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr0 --> v_addr25] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
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
void ref_p_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr1(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr1 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr1 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr1 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr1 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr1(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr1 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr1 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr1 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr1 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr1 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_v_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr2(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr2 --> v_addr0] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr2(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr2 --> v_addr2] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
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
                    if ( calAddrv_addr12( t, i) == calAddrv_addr2(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr2 --> v_addr12] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr2 --> v_addr14] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr2 --> v_addr16] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr2 --> v_addr23] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr2 --> v_addr24] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr2 --> v_addr25] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
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
void ref_q_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr3(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr3 --> q_addr3] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr3 --> q_addr9] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr3 --> q_addr11] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr3 --> q_addr15] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr3(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr3 --> q_addr20] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr3 --> q_addr26] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr3 --> q_addr28] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr3 --> q_addr32] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_v_addr12() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr12(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr12 --> v_addr0] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr12(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr12 --> v_addr2] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
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
                    if ( calAddrv_addr12( t, i) == calAddrv_addr12(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr12 --> v_addr12] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr12 --> v_addr14] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr12 --> v_addr16] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr12 --> v_addr23] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr12 --> v_addr24] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr12 --> v_addr25] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
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
void ref_p_addr13() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        if (j_Start % -1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr13 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr13 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr13 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr13 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr13 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr13 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr13 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr13 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr13 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_v_addr14() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        if (j_Start % -1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr14 --> v_addr0] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr14 --> v_addr2] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
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
                    if ( calAddrv_addr12( t, i) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr14 --> v_addr12] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr14 --> v_addr14] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr14 --> v_addr16] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr14 --> v_addr23] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr14 --> v_addr24] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr14 --> v_addr25] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
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
void ref_q_addr15() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        if (j_Start % -1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr15 --> q_addr3] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr15 --> q_addr9] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr15 --> q_addr11] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr15 --> q_addr15] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr15 --> q_addr20] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr15 --> q_addr26] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr15 --> q_addr28] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr15 --> q_addr32] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_p_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr5 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr5 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr5 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr5 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr5 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr5 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr5 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr5 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr5 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_u_addr17() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr17 --> u_addr6] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr17 --> u_addr7] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr17 --> u_addr8] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr17(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr17 --> u_addr17] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr17(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr17 --> u_addr19] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
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
                    if ( calAddru_addr29( t, i) == calAddru_addr17(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr17 --> u_addr29] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr17 --> u_addr31] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr17 --> u_addr33] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_u_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr6 --> u_addr6] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr6 --> u_addr7] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr6 --> u_addr8] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr6 --> u_addr17] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr6 --> u_addr19] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
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
                    if ( calAddru_addr29( t, i) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr6 --> u_addr29] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr6 --> u_addr31] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr6 --> u_addr33] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_q_addr20() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr20(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr20 --> q_addr3] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr20 --> q_addr9] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr20 --> q_addr11] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr20 --> q_addr15] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr20(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr20 --> q_addr20] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr20 --> q_addr26] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr20 --> q_addr28] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr20 --> q_addr32] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_p_addr21() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr21 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr21 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr21 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr21 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr21 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr21 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr21 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr21 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr21 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr21 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_p_addr22() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr22 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr22 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr22 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr22 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr22 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr22 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr22 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr22 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr22 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr22 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_v_addr23() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr23 --> v_addr0] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr23 --> v_addr2] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
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
                    if ( calAddrv_addr12( t, i) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr23 --> v_addr12] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr23 --> v_addr14] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr23 --> v_addr16] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr23 --> v_addr23] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr23 --> v_addr24] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr23 --> v_addr25] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
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
void ref_q_addr9() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr9 --> q_addr3] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr9 --> q_addr9] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr9 --> q_addr11] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr9 --> q_addr15] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr9 --> q_addr20] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr9 --> q_addr26] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr9 --> q_addr28] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr9 --> q_addr32] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_p_addr10() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr10 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr10 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr10 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr10 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr10 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr10 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr10 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr10 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr10 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_q_addr11() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr11 --> q_addr3] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr11 --> q_addr9] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr11 --> q_addr11] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr11 --> q_addr15] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr11 --> q_addr20] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr11 --> q_addr26] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr11 --> q_addr28] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr11 --> q_addr32] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_u_addr31() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        if (j_Start % -1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr31 --> u_addr6] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr31 --> u_addr7] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr31 --> u_addr8] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr31 --> u_addr17] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr31 --> u_addr19] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
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
                    if ( calAddru_addr29( t, i) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr31 --> u_addr29] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr31 --> u_addr31] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr31 --> u_addr33] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_q_addr32() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        if (j_Start % -1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr32 --> q_addr3] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr32 --> q_addr9] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr32 --> q_addr11] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr32 --> q_addr15] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr32 --> q_addr20] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr32 --> q_addr26] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr32 --> q_addr28] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr32 --> q_addr32] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_u_addr33() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        if (j_Start % -1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr33 --> u_addr6] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr33 --> u_addr7] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr33 --> u_addr8] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr33 --> u_addr17] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr33 --> u_addr19] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
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
                    if ( calAddru_addr29( t, i) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr33 --> u_addr29] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr33 --> u_addr31] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr33 --> u_addr33] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_v_addr16() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        if (j_Start % -1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr16 --> v_addr0] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr16 --> v_addr2] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
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
                    if ( calAddrv_addr12( t, i) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr16 --> v_addr12] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr16 --> v_addr14] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr16 --> v_addr16] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr16 --> v_addr23] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr16 --> v_addr24] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr16 --> v_addr25] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
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
void ref_p_addr18() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr18(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr18 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr18 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr18 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr18 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr18 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr18(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr18 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr18 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr18 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr18 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr18 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_u_addr19() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr19 --> u_addr6] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr19 --> u_addr7] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr19 --> u_addr8] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr19(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr19 --> u_addr17] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr19(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr19 --> u_addr19] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
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
                    if ( calAddru_addr29( t, i) == calAddru_addr19(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr19 --> u_addr29] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr19 --> u_addr31] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr19 --> u_addr33] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_v_addr24() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr24 --> v_addr0] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr24 --> v_addr2] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
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
                    if ( calAddrv_addr12( t, i) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr24 --> v_addr12] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr24 --> v_addr14] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr24 --> v_addr16] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr24 --> v_addr23] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr24 --> v_addr24] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr24 --> v_addr25] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
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
void ref_v_addr25() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr25 --> v_addr0] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr25 --> v_addr2] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
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
                    if ( calAddrv_addr12( t, i) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[v_addr25 --> v_addr12] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr25 --> v_addr14] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr25 --> v_addr16] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr25 --> v_addr23] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr25 --> v_addr24] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[v_addr25 --> v_addr25] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
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
void ref_q_addr26() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr26 --> q_addr3] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr26 --> q_addr9] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr26 --> q_addr11] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr26 --> q_addr15] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr26 --> q_addr20] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr26 --> q_addr26] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr26 --> q_addr28] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr26 --> q_addr32] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_p_addr27() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr27 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr27 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr27 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr27 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr27 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr27 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr27 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr27 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr27 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr27 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_q_addr28() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        if (j_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr28 --> q_addr3] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr28 --> q_addr9] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr28 --> q_addr11] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr28 --> q_addr15] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[q_addr28 --> q_addr20] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr28 --> q_addr26] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr28 --> q_addr28] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[q_addr28 --> q_addr32] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_u_addr29() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 102;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr29 --> u_addr6] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr29 --> u_addr7] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr29 --> u_addr8] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr29(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr29 --> u_addr17] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr29(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr29 --> u_addr19] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
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
                    if ( calAddru_addr29( t, i) == calAddru_addr29(t_Start, i_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1, i_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[u_addr29 --> u_addr29] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr29 --> u_addr31] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1, i_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[u_addr29 --> u_addr33] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
void ref_p_addr30() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10444;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if (t_Start % 1 != 0) goto SAMPLE; 
        if (t_Start + THREAD_NUM * CHUNK_SIZE > 10) { goto SAMPLE; }
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        if (j_Start % -1 != 0) goto SAMPLE; 
        string idx_string =  to_string(t_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t=(t + 1)) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr30 --> p_addr1] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr30 --> p_addr4] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr30 --> p_addr5] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr30 --> p_addr10] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr30 --> p_addr13] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i=(i + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1, i);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[p_addr30 --> p_addr18] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j=(j + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr30 --> p_addr21] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr30 --> p_addr22] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr30 --> p_addr27] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j=(j + -1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, i_Start, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, i, j);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((t_Start -1)) - 1) * CHUNK_SIZE * THREAD_NUM * 16721964 + getChunkID((t - 1)) * CHUNK_SIZE * THREAD_NUM * 16721964;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((t_Start -1), (t - 1), cnt, 16721964, 16721964, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[p_addr30 --> p_addr30] [" << parallel_rt << "] (" << t_Start<< ", " << i_Start<< ", " << j_Start<< ") --> (" << t<< ", " << i<< ", " << j<< ") " << endl;
#endif
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
int main() {
#ifdef PAPI_TIMER
    // Get starting timepoint
    auto start = high_resolution_clock::now();
#endif
    ref_p_addr4();
    ref_u_addr7();
    ref_u_addr8();
    ref_v_addr0();
    ref_p_addr1();
    ref_v_addr2();
    ref_q_addr3();
    ref_v_addr12();
    ref_p_addr13();
    ref_v_addr14();
    ref_q_addr15();
    ref_p_addr5();
    ref_u_addr17();
    ref_u_addr6();
    ref_q_addr20();
    ref_p_addr21();
    ref_p_addr22();
    ref_v_addr23();
    ref_q_addr9();
    ref_p_addr10();
    ref_q_addr11();
    ref_u_addr31();
    ref_q_addr32();
    ref_u_addr33();
    ref_v_addr16();
    ref_p_addr18();
    ref_u_addr19();
    ref_v_addr24();
    ref_v_addr25();
    ref_q_addr26();
    ref_p_addr27();
    ref_q_addr28();
    ref_u_addr29();
    ref_p_addr30();
    rtDump();
    RTtoMR_AET();
    dumpMR();
#ifdef PAPI_TIMER
    // Get ending timepoint
    auto stop = high_resolution_clock::now(); 
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
     cout << "Time taken by SPS:  " << duration.count() << endl; 
#endif
    return 0;
}
 /* Analyze function: adi */ 
