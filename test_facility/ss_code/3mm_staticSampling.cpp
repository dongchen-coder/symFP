
 /* Start to analysis array index
Array index info: Total number of references: 15
E.addr ((i * 128) + j)
E.addr ((i * 128) + j)
E.addr ((i * 128) + j)
C.addr ((i * 128) + k)
D.addr ((k * 128) + j)
A.addr ((i * 128) + k)
B.addr ((k * 128) + j)
G.addr ((i * 128) + j)
E.addr ((i * 128) + k)
F.addr ((k * 128) + j)
F.addr ((i * 128) + j)
F.addr ((i * 128) + j)
F.addr ((i * 128) + j)
G.addr ((i * 128) + j)
G.addr ((i * 128) + j)

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
--Loop Bound: (0, 128)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 128)
----Loop inc: (j + 1)
----Loop predicate: <
------array access E.addr ((i * 128) + j)
------k
------Loop Bound: (0, 128)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 128) + k)
--------array access B.addr ((k * 128) + j)
--------array access E.addr ((i * 128) + j)
--------array access E.addr ((i * 128) + j)
--i
--Loop Bound: (0, 128)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 128)
----Loop inc: (j + 1)
----Loop predicate: <
------array access F.addr ((i * 128) + j)
------k
------Loop Bound: (0, 128)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access C.addr ((i * 128) + k)
--------array access D.addr ((k * 128) + j)
--------array access F.addr ((i * 128) + j)
--------array access F.addr ((i * 128) + j)
--i
--Loop Bound: (0, 128)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 128)
----Loop inc: (j + 1)
----Loop predicate: <
------array access G.addr ((i * 128) + j)
------k
------Loop Bound: (0, 128)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access E.addr ((i * 128) + k)
--------array access F.addr ((k * 128) + j)
--------array access G.addr ((i * 128) + j)
--------array access G.addr ((i * 128) + j)

Finish analysis loops */ 
/* # of Out-most Loops: 3 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 128
------Sample number: 16384
--------Sample number: 2097152
----Sample number: 128
------Sample number: 16384
--------Sample number: 2097152
----Sample number: 128
------Sample number: 16384
--------Sample number: 2097152
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* E_addr3	8404992 */
/* E_addr4	8404992 */
/* E_addr0	8404992 */
/* C_addr6	8404992 */
/* D_addr7	8404992 */
/* A_addr1	8404992 */
/* B_addr2	8404992 */
/* G_addr10	8404992 */
/* E_addr11	8404992 */
/* F_addr12	8404992 */
/* F_addr5	8404992 */
/* F_addr8	8404992 */
/* F_addr9	8404992 */
/* G_addr13	8404992 */
/* G_addr14	8404992 */
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
#include <cmath>
#ifdef PAPI_TIMER
#  include <chrono>
#endif
using namespace std;
#ifdef PAPI_TIMER
using namespace  chrono;
#endif
 map<uint64_t, double> RT;
 map<uint64_t, double> MR;
int getChunkID(uint64_t i) {
    return floor(i / (CHUNK_SIZE * THREAD_NUM));
}
int getThreadID(uint64_t i) {
    return i / CHUNK_SIZE - floor(i / (CHUNK_SIZE * THREAD_NUM))*THREAD_NUM ;
}
int getThreadLocalPos(uint64_t i) {
    return i % CHUNK_SIZE;
}
uint64_t parallel_predict(uint64_t i_src, uint64_t i_sink, uint64_t rt, uint64_t lsrc, uint64_t lsink, bool is_normal_ref) {
    uint64_t parallel_rt = rt;
    if (!is_normal_ref && getThreadID(i_src) < THREAD_NUM-1) {
#ifdef DEBUG
            cout << "Neighboring Effect" << endl;
#endif
        return 1;
    }
    int tsrc = getThreadID(i_src);
    int tsink = getThreadID(i_sink);
    int dT = tsink - tsrc;
    /* intra chunk reuse */
    if (getChunkID(i_src) == getChunkID(i_sink)) {
        /* same thread -- scaling effect */
        if (dT == 0) {
#ifdef DEBUG
            cout << "Scaling Effect" << endl;
#endif
            parallel_rt = (rt * THREAD_NUM - 1);
        }
        else if (getThreadLocalPos(i_src) <= getThreadLocalPos(i_sink)) { // src-sink order
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
            parallel_rt = CHUNK_SIZE * lsrc * THREAD_NUM * dT - (rt * THREAD_NUM) - abs(dT);
        }
    } else { // inter chunk reuse 
#ifdef DEBUG
            cout << "Inter Chunk Reuse" << endl;
#endif
            parallel_rt = rt * THREAD_NUM - CHUNK_SIZE * THREAD_NUM * (lsrc*(THREAD_NUM - tsrc) + lsink * tsink) + CHUNK_SIZE * THREAD_NUM * lsink + dT;
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
/* Array E_addr	i j */ 
/* E_addr ((i * 128) + j) 0 */
int calAddrE_addr0( int i, int j) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* Array A_addr	i k */ 
/* A_addr ((i * 128) + k) 1 */
int calAddrA_addr1( int i, int j, int k) {
    int result = (((i * 128) + k)) * 8 / 64;
    return result;
}
/* Array B_addr	k j */ 
/* B_addr ((k * 128) + j) 2 */
int calAddrB_addr2( int i, int j, int k) {
    int result = (((k * 128) + j)) * 8 / 64;
    return result;
}
/* Array E_addr	i j */ 
/* E_addr ((i * 128) + j) 3 */
int calAddrE_addr3( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* Array E_addr	i j */ 
/* E_addr ((i * 128) + j) 4 */
int calAddrE_addr4( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* Array F_addr	i j */ 
/* F_addr ((i * 128) + j) 5 */
int calAddrF_addr5( int i, int j) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* Array C_addr	i k */ 
/* C_addr ((i * 128) + k) 6 */
int calAddrC_addr6( int i, int j, int k) {
    int result = (((i * 128) + k)) * 8 / 64;
    return result;
}
/* Array D_addr	k j */ 
/* D_addr ((k * 128) + j) 7 */
int calAddrD_addr7( int i, int j, int k) {
    int result = (((k * 128) + j)) * 8 / 64;
    return result;
}
/* Array F_addr	i j */ 
/* F_addr ((i * 128) + j) 8 */
int calAddrF_addr8( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* Array F_addr	i j */ 
/* F_addr ((i * 128) + j) 9 */
int calAddrF_addr9( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* Array G_addr	i j */ 
/* G_addr ((i * 128) + j) 10 */
int calAddrG_addr10( int i, int j) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* Array E_addr	i k */ 
/* E_addr ((i * 128) + k) 11 */
int calAddrE_addr11( int i, int j, int k) {
    int result = (((i * 128) + k)) * 8 / 64;
    return result;
}
/* Array F_addr	k j */ 
/* F_addr ((k * 128) + j) 12 */
int calAddrF_addr12( int i, int j, int k) {
    int result = (((k * 128) + j)) * 8 / 64;
    return result;
}
/* Array G_addr	i j */ 
/* G_addr ((i * 128) + j) 13 */
int calAddrG_addr13( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* Array G_addr	i j */ 
/* G_addr ((i * 128) + j) 14 */
int calAddrG_addr14( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
void ref_E_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 128; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 128; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrE_addr0( i, j) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                {
                int kLB2 = 0;
                if ( iLB0 == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr4( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 128; k=(k + 1)) {
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
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
void ref_E_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 128; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 128; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrE_addr0( i, j) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                {
                int kLB2 = 0;
                if ( iLB0 == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr4( i, j, k) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 128; k=(k + 1)) {
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
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
void ref_E_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16384;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 128; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 128; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrE_addr0( i, j) == calAddrE_addr0(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr4( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 128; k=(k + 1)) {
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
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
void ref_C_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( iLB3 == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                if ( iLB3 == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 128; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrC_addr6( i, j, k) == calAddrC_addr6(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        int iLB6 = 0;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
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
void ref_D_addr7() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( iLB3 == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                if ( iLB3 == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr7( i, j, k) == calAddrD_addr7(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, false);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        int iLB6 = 0;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
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
void ref_A_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 128; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( iLB0 == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 128; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 128; k=(k + 1)) {
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
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
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
void ref_B_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 128; i=(i + 1)) {
            {
            int jLB1 = 0;
            if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( iLB0 == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( i, j, k) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, false);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 128; k=(k + 1)) {
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
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
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
void ref_G_addr10() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16384;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( iLB6 == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrG_addr10( i, j) == calAddrG_addr10(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr13( i, j, k) == calAddrG_addr10(i_Start, j_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr14( i, j, k) == calAddrG_addr10(i_Start, j_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
void ref_E_addr11() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( iLB6 == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                if ( iLB6 == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr11(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
void ref_F_addr12() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( iLB6 == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                if ( iLB6 == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr12(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, false);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
void ref_F_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16384;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( iLB3 == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrF_addr5( i, j) == calAddrF_addr5(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr8( i, j, k) == calAddrF_addr5(i_Start, j_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr9( i, j, k) == calAddrF_addr5(i_Start, j_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        int iLB6 = 0;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr5(i_Start, j_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
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
void ref_F_addr8() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( iLB3 == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrF_addr5( i, j) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                if ( iLB3 == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr8( i, j, k) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr9( i, j, k) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        int iLB6 = 0;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
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
void ref_F_addr9() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 128; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( iLB3 == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 128; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrF_addr5( i, j) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                if ( iLB3 == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr8( i, j, k) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr9( i, j, k) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
        int iLB6 = 0;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
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
void ref_G_addr13() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( iLB6 == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrG_addr10( i, j) == calAddrG_addr13(i_Start, j_Start, k_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                {
                int kLB8 = 0;
                if ( iLB6 == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr13( i, j, k) == calAddrG_addr13(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr14( i, j, k) == calAddrG_addr13(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
void ref_G_addr14() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097152;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string =  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 128; i=(i + 1)) {
            {
            int jLB7 = 0;
            if ( iLB6 == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 128; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrG_addr10( i, j) == calAddrG_addr14(i_Start, j_Start, k_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                {
                int kLB8 = 0;
                if ( iLB6 == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 128; k=(k + 1)) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr13( i, j, k) == calAddrG_addr14(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
#endif
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr14( i, j, k) == calAddrG_addr14(i_Start, j_Start, k_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8404992, 8404992, true);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                if (parallel_rt != 1 && parallel_rt != 13) {
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ", " << k_Start<< ") -- (" << i<< ", " << j<< ", " << k<< ") " << endl;
                                }
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
int main() {
#ifdef PAPI_TIMER
    // Get starting timepoint
    auto start = high_resolution_clock::now();
#endif
    ref_E_addr3();
    ref_E_addr4();
    ref_E_addr0();
    ref_C_addr6();
    ref_D_addr7();
    ref_A_addr1();
    ref_B_addr2();
    ref_G_addr10();
    ref_E_addr11();
    ref_F_addr12();
    ref_F_addr5();
    ref_F_addr8();
    ref_F_addr9();
    ref_G_addr13();
    ref_G_addr14();
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
