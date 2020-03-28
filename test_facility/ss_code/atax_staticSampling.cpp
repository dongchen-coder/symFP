
 /* Start to analysis array index
Array index info: Total number of references: 8
tmp.addr i
tmp.addr i
y.addr j
A.addr ((i * 1024) + j)
x.addr j
tmp.addr i
y.addr j
A.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %nx
i32 %ny
double* %A
double* %x
double* %y
double* %tmp

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
------array access tmp.addr i
------array access A.addr ((i * 1024) + j)
------array access x.addr j
------array access tmp.addr i
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access y.addr j
------array access A.addr ((i * 1024) + j)
------array access tmp.addr i
------array access y.addr j

Finish analysis loops */ 
/* # of Out-most Loops: 1 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 1024
------Sample number: 1048576
------Sample number: 1048576
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* tmp_addr0	8388608 */
/* tmp_addr6	8388608 */
/* y_addr7	8388608 */
/* A_addr1	8388608 */
/* x_addr2	8388608 */
/* tmp_addr3	8388608 */
/* y_addr4	8388608 */
/* A_addr5	8388608 */
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
/* Array tmp_addr	i */ 
/* tmp_addr i 0 */
int calAddrtmp_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array A_addr	i j */ 
/* A_addr ((i * 1024) + j) 1 */
int calAddrA_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array x_addr	j */ 
/* x_addr j 2 */
int calAddrx_addr2( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* Array tmp_addr	i */ 
/* tmp_addr i 3 */
int calAddrtmp_addr3( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array y_addr	j */ 
/* y_addr j 4 */
int calAddry_addr4( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* Array A_addr	i j */ 
/* A_addr ((i * 1024) + j) 5 */
int calAddrA_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* Array tmp_addr	i */ 
/* tmp_addr i 6 */
int calAddrtmp_addr6( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array y_addr	j */ 
/* y_addr j 7 */
int calAddry_addr7( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
void ref_tmp_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1048576;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
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
            if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr0( i, j) == calAddrtmp_addr0(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
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
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr3( i, j) == calAddrtmp_addr0(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
            }
            }
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr6( i, j) == calAddrtmp_addr0(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
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
void ref_tmp_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1048576;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
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
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr0( i, j) == calAddrtmp_addr6(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr3( i, j) == calAddrtmp_addr6(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
            }
            }
            {
            int jLB2 = 0;
            if ( iLB0 == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr6( i, j) == calAddrtmp_addr6(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
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
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_y_addr7() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1048576;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
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
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB2 = 0;
            if ( iLB0 == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr4( i, j) == calAddry_addr7(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr7( i, j) == calAddry_addr7(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
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
    for ( int s = 0; s < 1048576;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
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
            if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( i, j) == calAddrA_addr1(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, true);
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
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr1(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_x_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1048576;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
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
            if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr2( i, j) == calAddrx_addr2(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
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
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 1024; j=(j + 1)) {
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
void ref_tmp_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1048576;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
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
            if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr0( i, j) == calAddrtmp_addr3(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr3( i, j) == calAddrtmp_addr3(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
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
            }
            }
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr6( i, j) == calAddrtmp_addr3(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
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
void ref_y_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1048576;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
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
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB2 = 0;
            if ( iLB0 == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr4( i, j) == calAddry_addr4(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
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
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr7( i, j) == calAddry_addr4(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, false);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
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
void ref_A_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1048576;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
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
            for ( int j = jLB1; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, true);
                            rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                            if (parallel_rt != 1 && parallel_rt != 13) {
                                cout << "[" << parallel_rt << "] (" << i_Start<< ", " << j_Start<< ") -- (" << i<< ", " << j<< ") " << endl;
                            }
#endif
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB2 = 0;
            if ( iLB0 == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                            uint64_t parallel_rt = parallel_predict((i_Start -0), i, cnt, 8388608, 8388608, true);
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
int main() {
#ifdef PAPI_TIMER
    // Get starting timepoint
    auto start = high_resolution_clock::now();
#endif
    ref_tmp_addr0();
    ref_tmp_addr6();
    ref_y_addr7();
    ref_A_addr1();
    ref_x_addr2();
    ref_tmp_addr3();
    ref_y_addr4();
    ref_A_addr5();
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
 /* Analyze function: atax */ 
