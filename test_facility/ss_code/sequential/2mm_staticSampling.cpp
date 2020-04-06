
 /* Start to analysis array index
Array index info: Total number of references: 11
tmp.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
C.addr ((k * 1024) + j)
D.addr ((i * 1024) + j)
D.addr ((i * 1024) + j)
B.addr ((k * 1024) + j)
tmp.addr ((i * 1024) + j)
tmp.addr ((i * 1024) + j)
D.addr ((i * 1024) + j)
D.addr ((i * 1024) + j)
tmp.addr ((i * 1024) + k)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %tmp
double* %A
double* %B
double* %C
double* %D
double %alpha
double %beta

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
------array access tmp.addr ((i * 1024) + j)
------k
------Loop Bound: (0, 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 1024) + k)
--------array access B.addr ((k * 1024) + j)
--------array access tmp.addr ((i * 1024) + j)
--------array access tmp.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access D.addr ((i * 1024) + j)
------array access D.addr ((i * 1024) + j)
------k
------Loop Bound: (0, 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access tmp.addr ((i * 1024) + k)
--------array access C.addr ((k * 1024) + j)
--------array access D.addr ((i * 1024) + j)
--------array access D.addr ((i * 1024) + j)

Finish analysis loops */ 
/* # of Out-most Loops: 2 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 102
------Sample number: 10485
--------Sample number: 1073741
----Sample number: 102
------Sample number: 10485
--------Sample number: 1073741
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* tmp_addr0	4195328 */
/* A_addr1	4195328 */
/* C_addr8	4196352 */
/* D_addr9	4196352 */
/* D_addr10	4196352 */
/* B_addr2	4195328 */
/* tmp_addr3	4195328 */
/* tmp_addr4	4195328 */
/* D_addr5	4196352 */
/* D_addr6	4196352 */
/* tmp_addr7	4196352 */
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
/* Array tmp_addr	i j */ 
/* tmp_addr ((i * 1024) + j) 0 */
int calAddrtmp_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array A_addr	i k */ 
/* A_addr ((i * 1024) + k) 1 */
int calAddrA_addr1( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* Array B_addr	k j */ 
/* B_addr ((k * 1024) + j) 2 */
int calAddrB_addr2( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
/* Array tmp_addr	i j */ 
/* tmp_addr ((i * 1024) + j) 3 */
int calAddrtmp_addr3( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array tmp_addr	i j */ 
/* tmp_addr ((i * 1024) + j) 4 */
int calAddrtmp_addr4( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array D_addr	i j */ 
/* D_addr ((i * 1024) + j) 5 */
int calAddrD_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array D_addr	i j */ 
/* D_addr ((i * 1024) + j) 6 */
int calAddrD_addr6( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array tmp_addr	i k */ 
/* tmp_addr ((i * 1024) + k) 7 */
int calAddrtmp_addr7( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* Array C_addr	k j */ 
/* C_addr ((k * 1024) + j) 8 */
int calAddrC_addr8( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
/* Array D_addr	i j */ 
/* D_addr ((i * 1024) + j) 9 */
int calAddrD_addr9( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array D_addr	i j */ 
/* D_addr ((i * 1024) + j) 10 */
int calAddrD_addr10( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_tmp_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10485;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                    if ( calAddrtmp_addr0( i, j) == calAddrtmp_addr0(i_Start, j_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr0, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr0, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[tmp_addr0 --> tmp_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
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
                        if ( calAddrtmp_addr3( i, j, k) == calAddrtmp_addr0(i_Start, j_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr0, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr3, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr0 --> tmp_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr4( i, j, k) == calAddrtmp_addr0(i_Start, j_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr0, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr4, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr0 --> tmp_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr7( i, j, k) == calAddrtmp_addr0(i_Start, j_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr0, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr7, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + (getChunkNum(1024) - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4196352, middle_accesses, true, false, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr0 --> tmp_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
void ref_A_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073741;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr1, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr1, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[A_addr1 --> A_addr1] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
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
EndSample:
        s++;
        }
}
void ref_C_addr8() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073741;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                        if ( calAddrC_addr8( i, j, k) == calAddrC_addr8(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrC_addr8, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrC_addr8, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[C_addr8 --> C_addr8] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
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
void ref_D_addr9() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073741;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                    if ( calAddrD_addr5( i, j) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr9, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr5, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[D_addr9 --> D_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr6( i, j) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr9, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr6, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[D_addr9 --> D_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
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
                        if ( calAddrD_addr9( i, j, k) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr9, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[D_addr9 --> D_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr10( i, j, k) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr9, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr10, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[D_addr9 --> D_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
void ref_D_addr10() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073741;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                    if ( calAddrD_addr5( i, j) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr10, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr5, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[D_addr10 --> D_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr6( i, j) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr10, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr6, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[D_addr10 --> D_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
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
                        if ( calAddrD_addr9( i, j, k) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr9, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[D_addr10 --> D_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr10( i, j, k) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr10, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr10, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[D_addr10 --> D_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
void ref_B_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073741;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                                /* is_normal_ref: false */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrB_addr2, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrB_addr2, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, false, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[B_addr2 --> B_addr2] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
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
EndSample:
        s++;
        }
}
void ref_tmp_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073741;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                    if ( calAddrtmp_addr0( i, j) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr3, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr0, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[tmp_addr3 --> tmp_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
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
                        if ( calAddrtmp_addr3( i, j, k) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr3, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr3 --> tmp_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr4( i, j, k) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr4, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr3 --> tmp_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr7( i, j, k) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr3, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr7, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + (getChunkNum(1024) - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4196352, middle_accesses, true, false, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr3 --> tmp_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
void ref_tmp_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073741;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                    if ( calAddrtmp_addr0( i, j) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr4, _1, j_Start, k_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr0, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[tmp_addr4 --> tmp_addr0] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
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
                        if ( calAddrtmp_addr3( i, j, k) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr3, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr4 --> tmp_addr3] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr4( i, j, k) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr4, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4195328;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4195328, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr4 --> tmp_addr4] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr7( i, j, k) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: false */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr4, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr7, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + (getChunkNum(1024) - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4195328 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4195328, 4196352, middle_accesses, true, false, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr4 --> tmp_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
void ref_D_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10485;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                    if ( calAddrD_addr5( i, j) == calAddrD_addr5(i_Start, j_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr5, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[D_addr5 --> D_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr6( i, j) == calAddrD_addr5(i_Start, j_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr5, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr6, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[D_addr5 --> D_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr9( i, j, k) == calAddrD_addr5(i_Start, j_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr5, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr9, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[D_addr5 --> D_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr10( i, j, k) == calAddrD_addr5(i_Start, j_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr5, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr10, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[D_addr5 --> D_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
void ref_D_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10485;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                    if ( calAddrD_addr5( i, j) == calAddrD_addr6(i_Start, j_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr5, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[D_addr6 --> D_addr5] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr6( i, j) == calAddrD_addr6(i_Start, j_Start)) {
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr6, _1, j_Start);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr6, _1, j);
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                            cout << " middle_access is " << middle_accesses << endl;
                            uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                            if (parallel_rt == 0) { goto EndSample; }
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                cout << "[D_addr6 --> D_addr6] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ") " << endl;
#endif
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
                        if ( calAddrD_addr9( i, j, k) == calAddrD_addr6(i_Start, j_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr6, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr9, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[D_addr6 --> D_addr9] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr10( i, j, k) == calAddrD_addr6(i_Start, j_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrD_addr6, _1, j_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrD_addr10, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[D_addr6 --> D_addr10] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
void ref_tmp_addr7() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073741;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if (i_Start + THREAD_NUM * CHUNK_SIZE > 1024) { goto SAMPLE; }
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
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 1024; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr7( i, j, k) == calAddrtmp_addr7(i_Start, j_Start, k_Start)) {
                                /* is_normal_ref: true */
                                /* is_in_same_loop: true */
                                /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                                function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrtmp_addr7, _1, j_Start, k_Start);
                                function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrtmp_addr7, _1, j, k);
                                /* compute the number of accesses between source and sink chunk */
                                uint64_t middle_accesses = 0 + ( - getChunkID((i_Start -0)) - 1) * CHUNK_SIZE * THREAD_NUM * 4196352 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 4196352;
                                cout << " middle_access is " << middle_accesses << endl;
                                uint64_t parallel_rt = parallel_predict((i_Start -0), (i - 0), cnt, 4196352, 4196352, middle_accesses, true, true, srcAddrCal, sinkAddrCal);
                                if (parallel_rt == 0) { goto EndSample; }
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[tmp_addr7 --> tmp_addr7] [" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") --> (" << i<< ", " << j<< ", " << k<< ") " << endl;
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
    ref_tmp_addr0();
    ref_A_addr1();
    ref_C_addr8();
    ref_D_addr9();
    ref_D_addr10();
    ref_B_addr2();
    ref_tmp_addr3();
    ref_tmp_addr4();
    ref_D_addr5();
    ref_D_addr6();
    ref_tmp_addr7();
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
 /* Analyze function: mm2 */ 
