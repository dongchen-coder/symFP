
 /* Start to analysis array index
Array index info
A.addr (((i * 1024) + j) - 1)
A.addr ((i * 1024) + j)
A.addr (((i * 1024) + 1) + j)
A.addr (((1 + i) * 1024) + j)
A.addr (((i - 1) * 1024) + j)
B.addr ((i * 1024) + j)
B.addr ((i * 1024) + j)
B.addr (((i * 1024) + j) - 1)
B.addr (((i * 1024) + 1) + j)
B.addr (((1 + i) * 1024) + j)
B.addr (((i - 1) * 1024) + j)
A.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %B

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (0, 10)
--Loop inc: (t + 1)
--Loop predicate: <
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (1, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access A.addr ((i * 1024) + j)
--------array access A.addr (((i * 1024) + j) - 1)
--------array access A.addr (((i * 1024) + 1) + j)
--------array access A.addr (((1 + i) * 1024) + j)
--------array access A.addr (((i - 1) * 1024) + j)
--------array access B.addr ((i * 1024) + j)
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (1, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access B.addr ((i * 1024) + j)
--------array access B.addr (((i * 1024) + j) - 1)
--------array access B.addr (((i * 1024) + 1) + j)
--------array access B.addr (((1 + i) * 1024) + j)
--------array access B.addr (((i - 1) * 1024) + j)
--------array access A.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 0
------Sample number: 0
--------Sample number: 0
------Sample number: 0
--------Sample number: 0
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
#include <cmath>
using namespace std;
std::map<uint64_t, double> RT;
std::map<uint64_t, double> MR;
void rtHistoCal( int rt) {
    if (RT.find(rt) == RT.end()) { 
        RT[rt] = 1;
    } else {
        RT[rt] += 1;
    }
    return;
}
void RTtoMR_AET() {
    std::map<uint64_t, double> P;
    double total_num_RT = 0;
    uint64_t max_RT = 0;
    for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
        total_num_RT += it->second;
        if (max_RT < it->first) {
            max_RT = it->first;
        }
    }
    double accumulate_num_RT = 0;
    for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
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
        cout << it->first << " " << it->second << "\n";
    }
    return;
}
void dumpMR() {
    cout << "miss ratio" << endl;
    std::map<uint64_t, double>::iterator it1 = MR.begin();
    std::map<uint64_t, double>::iterator it2 = MR.begin();
    while(it1 != MR.end()) {
        while(1) {
            std::map<uint64_t, double>::iterator it3 = it2;
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
        cout << it1->first << " " << it1->second << endl;
        if (it1 != it2) {
            cout << it2->first << " " << it2->second << endl;
        }
        it1 = ++it2;
        it2 = it1;
    }
    return;
}
int calAddrA_addr0( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr1( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrA_addr2( int t, int i, int j) {
    int result = ((((i * 1024) + 1) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr3( int t, int i, int j) {
    int result = ((((1 + i) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr4( int t, int i, int j) {
    int result = ((((i - 1) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrB_addr0( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrB_addr1( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrB_addr2( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrB_addr3( int t, int i, int j) {
    int result = ((((i * 1024) + 1) + j)) * 8 / 64;
    return result;
}
int calAddrB_addr4( int t, int i, int j) {
    int result = ((((1 + i) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrB_addr5( int t, int i, int j) {
    int result = ((((i - 1) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr5( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_t_Start_A_addr0 = -1;
    uint64_t prev_t_End_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_t_Start_A_addr1 = -1;
    uint64_t prev_t_End_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_t_Start_A_addr2 = -1;
    uint64_t prev_t_End_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_t_Start_A_addr3 = -1;
    uint64_t prev_t_End_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_t_Start_A_addr4 = -1;
    uint64_t prev_t_End_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( t_Start - prev_t_Start_A_addr0 + prev_t_End_A_addr0, i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( t_Start - prev_t_Start_A_addr1 + prev_t_End_A_addr1, i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( t_Start - prev_t_Start_A_addr2 + prev_t_End_A_addr2, i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( t_Start - prev_t_Start_A_addr3 + prev_t_End_A_addr3, i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( t_Start - prev_t_Start_A_addr4 + prev_t_End_A_addr4, i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_t_Start_A_addr0 = t_Start;
                            prev_t_End_A_addr0 = t;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_t_Start_A_addr1 = t_Start;
                            prev_t_End_A_addr1 = t;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr2 = cnt;
                            prev_t_Start_A_addr2 = t_Start;
                            prev_t_End_A_addr2 = t;
                            prev_i_Start_A_addr2 = i_Start;
                            prev_i_End_A_addr2 = i;
                            prev_j_Start_A_addr2 = j_Start;
                            prev_j_End_A_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_t_Start_A_addr3 = t_Start;
                            prev_t_End_A_addr3 = t;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr4 = cnt;
                            prev_t_Start_A_addr4 = t_Start;
                            prev_t_End_A_addr4 = t;
                            prev_i_Start_A_addr4 = i_Start;
                            prev_i_End_A_addr4 = i;
                            prev_j_Start_A_addr4 = j_Start;
                            prev_j_End_A_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB3 = 1;
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
void ref_A_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_t_Start_A_addr0 = -1;
    uint64_t prev_t_End_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_t_Start_A_addr1 = -1;
    uint64_t prev_t_End_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_t_Start_A_addr2 = -1;
    uint64_t prev_t_End_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_t_Start_A_addr3 = -1;
    uint64_t prev_t_End_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_t_Start_A_addr4 = -1;
    uint64_t prev_t_End_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( t_Start - prev_t_Start_A_addr0 + prev_t_End_A_addr0, i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( t_Start - prev_t_Start_A_addr1 + prev_t_End_A_addr1, i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( t_Start - prev_t_Start_A_addr2 + prev_t_End_A_addr2, i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( t_Start - prev_t_Start_A_addr3 + prev_t_End_A_addr3, i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( t_Start - prev_t_Start_A_addr4 + prev_t_End_A_addr4, i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_t_Start_A_addr0 = t_Start;
                            prev_t_End_A_addr0 = t;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_t_Start_A_addr1 = t_Start;
                            prev_t_End_A_addr1 = t;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr2 = cnt;
                            prev_t_Start_A_addr2 = t_Start;
                            prev_t_End_A_addr2 = t;
                            prev_i_Start_A_addr2 = i_Start;
                            prev_i_End_A_addr2 = i;
                            prev_j_Start_A_addr2 = j_Start;
                            prev_j_End_A_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_t_Start_A_addr3 = t_Start;
                            prev_t_End_A_addr3 = t;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr4 = cnt;
                            prev_t_Start_A_addr4 = t_Start;
                            prev_t_End_A_addr4 = t;
                            prev_i_Start_A_addr4 = i_Start;
                            prev_i_End_A_addr4 = i;
                            prev_j_Start_A_addr4 = j_Start;
                            prev_j_End_A_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB3 = 1;
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
void ref_A_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_t_Start_A_addr0 = -1;
    uint64_t prev_t_End_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_t_Start_A_addr1 = -1;
    uint64_t prev_t_End_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_t_Start_A_addr2 = -1;
    uint64_t prev_t_End_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_t_Start_A_addr3 = -1;
    uint64_t prev_t_End_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_t_Start_A_addr4 = -1;
    uint64_t prev_t_End_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( t_Start - prev_t_Start_A_addr0 + prev_t_End_A_addr0, i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( t_Start - prev_t_Start_A_addr1 + prev_t_End_A_addr1, i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( t_Start - prev_t_Start_A_addr2 + prev_t_End_A_addr2, i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( t_Start - prev_t_Start_A_addr3 + prev_t_End_A_addr3, i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( t_Start - prev_t_Start_A_addr4 + prev_t_End_A_addr4, i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_t_Start_A_addr0 = t_Start;
                            prev_t_End_A_addr0 = t;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_t_Start_A_addr1 = t_Start;
                            prev_t_End_A_addr1 = t;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr2 = cnt;
                            prev_t_Start_A_addr2 = t_Start;
                            prev_t_End_A_addr2 = t;
                            prev_i_Start_A_addr2 = i_Start;
                            prev_i_End_A_addr2 = i;
                            prev_j_Start_A_addr2 = j_Start;
                            prev_j_End_A_addr2 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_t_Start_A_addr3 = t_Start;
                            prev_t_End_A_addr3 = t;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr4 = cnt;
                            prev_t_Start_A_addr4 = t_Start;
                            prev_t_End_A_addr4 = t;
                            prev_i_Start_A_addr4 = i_Start;
                            prev_i_End_A_addr4 = i;
                            prev_j_Start_A_addr4 = j_Start;
                            prev_j_End_A_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB3 = 1;
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
void ref_A_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_t_Start_A_addr0 = -1;
    uint64_t prev_t_End_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_t_Start_A_addr1 = -1;
    uint64_t prev_t_End_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_t_Start_A_addr2 = -1;
    uint64_t prev_t_End_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_t_Start_A_addr3 = -1;
    uint64_t prev_t_End_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_t_Start_A_addr4 = -1;
    uint64_t prev_t_End_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( t_Start - prev_t_Start_A_addr0 + prev_t_End_A_addr0, i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( t_Start - prev_t_Start_A_addr1 + prev_t_End_A_addr1, i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( t_Start - prev_t_Start_A_addr2 + prev_t_End_A_addr2, i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( t_Start - prev_t_Start_A_addr3 + prev_t_End_A_addr3, i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( t_Start - prev_t_Start_A_addr4 + prev_t_End_A_addr4, i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_t_Start_A_addr0 = t_Start;
                            prev_t_End_A_addr0 = t;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_t_Start_A_addr1 = t_Start;
                            prev_t_End_A_addr1 = t;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr2 = cnt;
                            prev_t_Start_A_addr2 = t_Start;
                            prev_t_End_A_addr2 = t;
                            prev_i_Start_A_addr2 = i_Start;
                            prev_i_End_A_addr2 = i;
                            prev_j_Start_A_addr2 = j_Start;
                            prev_j_End_A_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_t_Start_A_addr3 = t_Start;
                            prev_t_End_A_addr3 = t;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr4 = cnt;
                            prev_t_Start_A_addr4 = t_Start;
                            prev_t_End_A_addr4 = t;
                            prev_i_Start_A_addr4 = i_Start;
                            prev_i_End_A_addr4 = i;
                            prev_j_Start_A_addr4 = j_Start;
                            prev_j_End_A_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB3 = 1;
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
void ref_A_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_t_Start_A_addr0 = -1;
    uint64_t prev_t_End_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_t_Start_A_addr1 = -1;
    uint64_t prev_t_End_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_t_Start_A_addr2 = -1;
    uint64_t prev_t_End_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_t_Start_A_addr3 = -1;
    uint64_t prev_t_End_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_t_Start_A_addr4 = -1;
    uint64_t prev_t_End_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( t_Start - prev_t_Start_A_addr0 + prev_t_End_A_addr0, i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( t_Start - prev_t_Start_A_addr1 + prev_t_End_A_addr1, i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( t_Start - prev_t_Start_A_addr2 + prev_t_End_A_addr2, i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( t_Start - prev_t_Start_A_addr3 + prev_t_End_A_addr3, i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( t_Start - prev_t_Start_A_addr4 + prev_t_End_A_addr4, i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_t_Start_A_addr0 = t_Start;
                            prev_t_End_A_addr0 = t;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_t_Start_A_addr1 = t_Start;
                            prev_t_End_A_addr1 = t;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr2 = cnt;
                            prev_t_Start_A_addr2 = t_Start;
                            prev_t_End_A_addr2 = t;
                            prev_i_Start_A_addr2 = i_Start;
                            prev_i_End_A_addr2 = i;
                            prev_j_Start_A_addr2 = j_Start;
                            prev_j_End_A_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_t_Start_A_addr3 = t_Start;
                            prev_t_End_A_addr3 = t;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr4 = cnt;
                            prev_t_Start_A_addr4 = t_Start;
                            prev_t_End_A_addr4 = t;
                            prev_i_Start_A_addr4 = i_Start;
                            prev_i_End_A_addr4 = i;
                            prev_j_Start_A_addr4 = j_Start;
                            prev_j_End_A_addr4 = j;
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
            int iLB3 = 1;
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
void ref_A_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr5 = -1;
    uint64_t prev_t_Start_A_addr5 = -1;
    uint64_t prev_t_End_A_addr5 = -1;
    uint64_t prev_i_Start_A_addr5 = -1;
    uint64_t prev_i_End_A_addr5 = -1;
    uint64_t prev_j_Start_A_addr5 = -1;
    uint64_t prev_j_End_A_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr5 != -1) {
            if ( calAddrA_addr5( t_Start - prev_t_Start_A_addr5 + prev_t_End_A_addr5, i_Start - prev_i_Start_A_addr5 + prev_i_End_A_addr5, j_Start - prev_j_Start_A_addr5 + prev_j_End_A_addr5) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB3 = 1;
            if ( t == t_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr5 = cnt;
                            prev_t_Start_A_addr5 = t_Start;
                            prev_t_End_A_addr5 = t;
                            prev_i_Start_A_addr5 = i_Start;
                            prev_i_End_A_addr5 = i;
                            prev_j_Start_A_addr5 = j_Start;
                            prev_j_End_A_addr5 = j;
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
void ref_B_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr0 = -1;
    uint64_t prev_t_Start_B_addr0 = -1;
    uint64_t prev_t_End_B_addr0 = -1;
    uint64_t prev_i_Start_B_addr0 = -1;
    uint64_t prev_i_End_B_addr0 = -1;
    uint64_t prev_j_Start_B_addr0 = -1;
    uint64_t prev_j_End_B_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr0 != -1) {
            if ( calAddrB_addr0( t_Start - prev_t_Start_B_addr0 + prev_t_End_B_addr0, i_Start - prev_i_Start_B_addr0 + prev_i_End_B_addr0, j_Start - prev_j_Start_B_addr0 + prev_j_End_B_addr0) == calAddrB_addr0(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( t, i, j) == calAddrB_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr0 = cnt;
                            prev_t_Start_B_addr0 = t_Start;
                            prev_t_End_B_addr0 = t;
                            prev_i_Start_B_addr0 = i_Start;
                            prev_i_End_B_addr0 = i;
                            prev_j_Start_B_addr0 = j_Start;
                            prev_j_End_B_addr0 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
            }
            }
            {
            int iLB3 = 1;
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( t, i, j) == calAddrB_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( t, i, j) == calAddrB_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr3( t, i, j) == calAddrB_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr4( t, i, j) == calAddrB_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr5( t, i, j) == calAddrB_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
void ref_B_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr1 = -1;
    uint64_t prev_t_Start_B_addr1 = -1;
    uint64_t prev_t_End_B_addr1 = -1;
    uint64_t prev_i_Start_B_addr1 = -1;
    uint64_t prev_i_End_B_addr1 = -1;
    uint64_t prev_j_Start_B_addr1 = -1;
    uint64_t prev_j_End_B_addr1 = -1;
    uint64_t prev_cnt_B_addr2 = -1;
    uint64_t prev_t_Start_B_addr2 = -1;
    uint64_t prev_t_End_B_addr2 = -1;
    uint64_t prev_i_Start_B_addr2 = -1;
    uint64_t prev_i_End_B_addr2 = -1;
    uint64_t prev_j_Start_B_addr2 = -1;
    uint64_t prev_j_End_B_addr2 = -1;
    uint64_t prev_cnt_B_addr3 = -1;
    uint64_t prev_t_Start_B_addr3 = -1;
    uint64_t prev_t_End_B_addr3 = -1;
    uint64_t prev_i_Start_B_addr3 = -1;
    uint64_t prev_i_End_B_addr3 = -1;
    uint64_t prev_j_Start_B_addr3 = -1;
    uint64_t prev_j_End_B_addr3 = -1;
    uint64_t prev_cnt_B_addr4 = -1;
    uint64_t prev_t_Start_B_addr4 = -1;
    uint64_t prev_t_End_B_addr4 = -1;
    uint64_t prev_i_Start_B_addr4 = -1;
    uint64_t prev_i_End_B_addr4 = -1;
    uint64_t prev_j_Start_B_addr4 = -1;
    uint64_t prev_j_End_B_addr4 = -1;
    uint64_t prev_cnt_B_addr5 = -1;
    uint64_t prev_t_Start_B_addr5 = -1;
    uint64_t prev_t_End_B_addr5 = -1;
    uint64_t prev_i_Start_B_addr5 = -1;
    uint64_t prev_i_End_B_addr5 = -1;
    uint64_t prev_j_Start_B_addr5 = -1;
    uint64_t prev_j_End_B_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr1 != -1) {
            if ( calAddrB_addr1( t_Start - prev_t_Start_B_addr1 + prev_t_End_B_addr1, i_Start - prev_i_Start_B_addr1 + prev_i_End_B_addr1, j_Start - prev_j_Start_B_addr1 + prev_j_End_B_addr1) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr2 != -1) {
            if ( calAddrB_addr2( t_Start - prev_t_Start_B_addr2 + prev_t_End_B_addr2, i_Start - prev_i_Start_B_addr2 + prev_i_End_B_addr2, j_Start - prev_j_Start_B_addr2 + prev_j_End_B_addr2) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr3 != -1) {
            if ( calAddrB_addr3( t_Start - prev_t_Start_B_addr3 + prev_t_End_B_addr3, i_Start - prev_i_Start_B_addr3 + prev_i_End_B_addr3, j_Start - prev_j_Start_B_addr3 + prev_j_End_B_addr3) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr4 != -1) {
            if ( calAddrB_addr4( t_Start - prev_t_Start_B_addr4 + prev_t_End_B_addr4, i_Start - prev_i_Start_B_addr4 + prev_i_End_B_addr4, j_Start - prev_j_Start_B_addr4 + prev_j_End_B_addr4) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr5 != -1) {
            if ( calAddrB_addr5( t_Start - prev_t_Start_B_addr5 + prev_t_End_B_addr5, i_Start - prev_i_Start_B_addr5 + prev_i_End_B_addr5, j_Start - prev_j_Start_B_addr5 + prev_j_End_B_addr5) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( t, i, j) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB3 = 1;
            if ( t == t_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( t, i, j) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr1 = cnt;
                            prev_t_Start_B_addr1 = t_Start;
                            prev_t_End_B_addr1 = t;
                            prev_i_Start_B_addr1 = i_Start;
                            prev_i_End_B_addr1 = i;
                            prev_j_Start_B_addr1 = j_Start;
                            prev_j_End_B_addr1 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( t, i, j) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr2 = cnt;
                            prev_t_Start_B_addr2 = t_Start;
                            prev_t_End_B_addr2 = t;
                            prev_i_Start_B_addr2 = i_Start;
                            prev_i_End_B_addr2 = i;
                            prev_j_Start_B_addr2 = j_Start;
                            prev_j_End_B_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr3( t, i, j) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr3 = cnt;
                            prev_t_Start_B_addr3 = t_Start;
                            prev_t_End_B_addr3 = t;
                            prev_i_Start_B_addr3 = i_Start;
                            prev_i_End_B_addr3 = i;
                            prev_j_Start_B_addr3 = j_Start;
                            prev_j_End_B_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr4( t, i, j) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr4 = cnt;
                            prev_t_Start_B_addr4 = t_Start;
                            prev_t_End_B_addr4 = t;
                            prev_i_Start_B_addr4 = i_Start;
                            prev_i_End_B_addr4 = i;
                            prev_j_Start_B_addr4 = j_Start;
                            prev_j_End_B_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr5( t, i, j) == calAddrB_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr5 = cnt;
                            prev_t_Start_B_addr5 = t_Start;
                            prev_t_End_B_addr5 = t;
                            prev_i_Start_B_addr5 = i_Start;
                            prev_i_End_B_addr5 = i;
                            prev_j_Start_B_addr5 = j_Start;
                            prev_j_End_B_addr5 = j;
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
void ref_B_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr1 = -1;
    uint64_t prev_t_Start_B_addr1 = -1;
    uint64_t prev_t_End_B_addr1 = -1;
    uint64_t prev_i_Start_B_addr1 = -1;
    uint64_t prev_i_End_B_addr1 = -1;
    uint64_t prev_j_Start_B_addr1 = -1;
    uint64_t prev_j_End_B_addr1 = -1;
    uint64_t prev_cnt_B_addr2 = -1;
    uint64_t prev_t_Start_B_addr2 = -1;
    uint64_t prev_t_End_B_addr2 = -1;
    uint64_t prev_i_Start_B_addr2 = -1;
    uint64_t prev_i_End_B_addr2 = -1;
    uint64_t prev_j_Start_B_addr2 = -1;
    uint64_t prev_j_End_B_addr2 = -1;
    uint64_t prev_cnt_B_addr3 = -1;
    uint64_t prev_t_Start_B_addr3 = -1;
    uint64_t prev_t_End_B_addr3 = -1;
    uint64_t prev_i_Start_B_addr3 = -1;
    uint64_t prev_i_End_B_addr3 = -1;
    uint64_t prev_j_Start_B_addr3 = -1;
    uint64_t prev_j_End_B_addr3 = -1;
    uint64_t prev_cnt_B_addr4 = -1;
    uint64_t prev_t_Start_B_addr4 = -1;
    uint64_t prev_t_End_B_addr4 = -1;
    uint64_t prev_i_Start_B_addr4 = -1;
    uint64_t prev_i_End_B_addr4 = -1;
    uint64_t prev_j_Start_B_addr4 = -1;
    uint64_t prev_j_End_B_addr4 = -1;
    uint64_t prev_cnt_B_addr5 = -1;
    uint64_t prev_t_Start_B_addr5 = -1;
    uint64_t prev_t_End_B_addr5 = -1;
    uint64_t prev_i_Start_B_addr5 = -1;
    uint64_t prev_i_End_B_addr5 = -1;
    uint64_t prev_j_Start_B_addr5 = -1;
    uint64_t prev_j_End_B_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr1 != -1) {
            if ( calAddrB_addr1( t_Start - prev_t_Start_B_addr1 + prev_t_End_B_addr1, i_Start - prev_i_Start_B_addr1 + prev_i_End_B_addr1, j_Start - prev_j_Start_B_addr1 + prev_j_End_B_addr1) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr2 != -1) {
            if ( calAddrB_addr2( t_Start - prev_t_Start_B_addr2 + prev_t_End_B_addr2, i_Start - prev_i_Start_B_addr2 + prev_i_End_B_addr2, j_Start - prev_j_Start_B_addr2 + prev_j_End_B_addr2) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr3 != -1) {
            if ( calAddrB_addr3( t_Start - prev_t_Start_B_addr3 + prev_t_End_B_addr3, i_Start - prev_i_Start_B_addr3 + prev_i_End_B_addr3, j_Start - prev_j_Start_B_addr3 + prev_j_End_B_addr3) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr4 != -1) {
            if ( calAddrB_addr4( t_Start - prev_t_Start_B_addr4 + prev_t_End_B_addr4, i_Start - prev_i_Start_B_addr4 + prev_i_End_B_addr4, j_Start - prev_j_Start_B_addr4 + prev_j_End_B_addr4) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr5 != -1) {
            if ( calAddrB_addr5( t_Start - prev_t_Start_B_addr5 + prev_t_End_B_addr5, i_Start - prev_i_Start_B_addr5 + prev_i_End_B_addr5, j_Start - prev_j_Start_B_addr5 + prev_j_End_B_addr5) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( t, i, j) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB3 = 1;
            if ( t == t_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( t, i, j) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr1 = cnt;
                            prev_t_Start_B_addr1 = t_Start;
                            prev_t_End_B_addr1 = t;
                            prev_i_Start_B_addr1 = i_Start;
                            prev_i_End_B_addr1 = i;
                            prev_j_Start_B_addr1 = j_Start;
                            prev_j_End_B_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( t, i, j) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr2 = cnt;
                            prev_t_Start_B_addr2 = t_Start;
                            prev_t_End_B_addr2 = t;
                            prev_i_Start_B_addr2 = i_Start;
                            prev_i_End_B_addr2 = i;
                            prev_j_Start_B_addr2 = j_Start;
                            prev_j_End_B_addr2 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr3( t, i, j) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr3 = cnt;
                            prev_t_Start_B_addr3 = t_Start;
                            prev_t_End_B_addr3 = t;
                            prev_i_Start_B_addr3 = i_Start;
                            prev_i_End_B_addr3 = i;
                            prev_j_Start_B_addr3 = j_Start;
                            prev_j_End_B_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr4( t, i, j) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr4 = cnt;
                            prev_t_Start_B_addr4 = t_Start;
                            prev_t_End_B_addr4 = t;
                            prev_i_Start_B_addr4 = i_Start;
                            prev_i_End_B_addr4 = i;
                            prev_j_Start_B_addr4 = j_Start;
                            prev_j_End_B_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr5( t, i, j) == calAddrB_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr5 = cnt;
                            prev_t_Start_B_addr5 = t_Start;
                            prev_t_End_B_addr5 = t;
                            prev_i_Start_B_addr5 = i_Start;
                            prev_i_End_B_addr5 = i;
                            prev_j_Start_B_addr5 = j_Start;
                            prev_j_End_B_addr5 = j;
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
void ref_B_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr1 = -1;
    uint64_t prev_t_Start_B_addr1 = -1;
    uint64_t prev_t_End_B_addr1 = -1;
    uint64_t prev_i_Start_B_addr1 = -1;
    uint64_t prev_i_End_B_addr1 = -1;
    uint64_t prev_j_Start_B_addr1 = -1;
    uint64_t prev_j_End_B_addr1 = -1;
    uint64_t prev_cnt_B_addr2 = -1;
    uint64_t prev_t_Start_B_addr2 = -1;
    uint64_t prev_t_End_B_addr2 = -1;
    uint64_t prev_i_Start_B_addr2 = -1;
    uint64_t prev_i_End_B_addr2 = -1;
    uint64_t prev_j_Start_B_addr2 = -1;
    uint64_t prev_j_End_B_addr2 = -1;
    uint64_t prev_cnt_B_addr3 = -1;
    uint64_t prev_t_Start_B_addr3 = -1;
    uint64_t prev_t_End_B_addr3 = -1;
    uint64_t prev_i_Start_B_addr3 = -1;
    uint64_t prev_i_End_B_addr3 = -1;
    uint64_t prev_j_Start_B_addr3 = -1;
    uint64_t prev_j_End_B_addr3 = -1;
    uint64_t prev_cnt_B_addr4 = -1;
    uint64_t prev_t_Start_B_addr4 = -1;
    uint64_t prev_t_End_B_addr4 = -1;
    uint64_t prev_i_Start_B_addr4 = -1;
    uint64_t prev_i_End_B_addr4 = -1;
    uint64_t prev_j_Start_B_addr4 = -1;
    uint64_t prev_j_End_B_addr4 = -1;
    uint64_t prev_cnt_B_addr5 = -1;
    uint64_t prev_t_Start_B_addr5 = -1;
    uint64_t prev_t_End_B_addr5 = -1;
    uint64_t prev_i_Start_B_addr5 = -1;
    uint64_t prev_i_End_B_addr5 = -1;
    uint64_t prev_j_Start_B_addr5 = -1;
    uint64_t prev_j_End_B_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr1 != -1) {
            if ( calAddrB_addr1( t_Start - prev_t_Start_B_addr1 + prev_t_End_B_addr1, i_Start - prev_i_Start_B_addr1 + prev_i_End_B_addr1, j_Start - prev_j_Start_B_addr1 + prev_j_End_B_addr1) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr2 != -1) {
            if ( calAddrB_addr2( t_Start - prev_t_Start_B_addr2 + prev_t_End_B_addr2, i_Start - prev_i_Start_B_addr2 + prev_i_End_B_addr2, j_Start - prev_j_Start_B_addr2 + prev_j_End_B_addr2) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr3 != -1) {
            if ( calAddrB_addr3( t_Start - prev_t_Start_B_addr3 + prev_t_End_B_addr3, i_Start - prev_i_Start_B_addr3 + prev_i_End_B_addr3, j_Start - prev_j_Start_B_addr3 + prev_j_End_B_addr3) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr4 != -1) {
            if ( calAddrB_addr4( t_Start - prev_t_Start_B_addr4 + prev_t_End_B_addr4, i_Start - prev_i_Start_B_addr4 + prev_i_End_B_addr4, j_Start - prev_j_Start_B_addr4 + prev_j_End_B_addr4) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr5 != -1) {
            if ( calAddrB_addr5( t_Start - prev_t_Start_B_addr5 + prev_t_End_B_addr5, i_Start - prev_i_Start_B_addr5 + prev_i_End_B_addr5, j_Start - prev_j_Start_B_addr5 + prev_j_End_B_addr5) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( t, i, j) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB3 = 1;
            if ( t == t_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( t, i, j) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr1 = cnt;
                            prev_t_Start_B_addr1 = t_Start;
                            prev_t_End_B_addr1 = t;
                            prev_i_Start_B_addr1 = i_Start;
                            prev_i_End_B_addr1 = i;
                            prev_j_Start_B_addr1 = j_Start;
                            prev_j_End_B_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( t, i, j) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr2 = cnt;
                            prev_t_Start_B_addr2 = t_Start;
                            prev_t_End_B_addr2 = t;
                            prev_i_Start_B_addr2 = i_Start;
                            prev_i_End_B_addr2 = i;
                            prev_j_Start_B_addr2 = j_Start;
                            prev_j_End_B_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr3( t, i, j) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr3 = cnt;
                            prev_t_Start_B_addr3 = t_Start;
                            prev_t_End_B_addr3 = t;
                            prev_i_Start_B_addr3 = i_Start;
                            prev_i_End_B_addr3 = i;
                            prev_j_Start_B_addr3 = j_Start;
                            prev_j_End_B_addr3 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr4( t, i, j) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr4 = cnt;
                            prev_t_Start_B_addr4 = t_Start;
                            prev_t_End_B_addr4 = t;
                            prev_i_Start_B_addr4 = i_Start;
                            prev_i_End_B_addr4 = i;
                            prev_j_Start_B_addr4 = j_Start;
                            prev_j_End_B_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr5( t, i, j) == calAddrB_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr5 = cnt;
                            prev_t_Start_B_addr5 = t_Start;
                            prev_t_End_B_addr5 = t;
                            prev_i_Start_B_addr5 = i_Start;
                            prev_i_End_B_addr5 = i;
                            prev_j_Start_B_addr5 = j_Start;
                            prev_j_End_B_addr5 = j;
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
void ref_B_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr1 = -1;
    uint64_t prev_t_Start_B_addr1 = -1;
    uint64_t prev_t_End_B_addr1 = -1;
    uint64_t prev_i_Start_B_addr1 = -1;
    uint64_t prev_i_End_B_addr1 = -1;
    uint64_t prev_j_Start_B_addr1 = -1;
    uint64_t prev_j_End_B_addr1 = -1;
    uint64_t prev_cnt_B_addr2 = -1;
    uint64_t prev_t_Start_B_addr2 = -1;
    uint64_t prev_t_End_B_addr2 = -1;
    uint64_t prev_i_Start_B_addr2 = -1;
    uint64_t prev_i_End_B_addr2 = -1;
    uint64_t prev_j_Start_B_addr2 = -1;
    uint64_t prev_j_End_B_addr2 = -1;
    uint64_t prev_cnt_B_addr3 = -1;
    uint64_t prev_t_Start_B_addr3 = -1;
    uint64_t prev_t_End_B_addr3 = -1;
    uint64_t prev_i_Start_B_addr3 = -1;
    uint64_t prev_i_End_B_addr3 = -1;
    uint64_t prev_j_Start_B_addr3 = -1;
    uint64_t prev_j_End_B_addr3 = -1;
    uint64_t prev_cnt_B_addr4 = -1;
    uint64_t prev_t_Start_B_addr4 = -1;
    uint64_t prev_t_End_B_addr4 = -1;
    uint64_t prev_i_Start_B_addr4 = -1;
    uint64_t prev_i_End_B_addr4 = -1;
    uint64_t prev_j_Start_B_addr4 = -1;
    uint64_t prev_j_End_B_addr4 = -1;
    uint64_t prev_cnt_B_addr5 = -1;
    uint64_t prev_t_Start_B_addr5 = -1;
    uint64_t prev_t_End_B_addr5 = -1;
    uint64_t prev_i_Start_B_addr5 = -1;
    uint64_t prev_i_End_B_addr5 = -1;
    uint64_t prev_j_Start_B_addr5 = -1;
    uint64_t prev_j_End_B_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr1 != -1) {
            if ( calAddrB_addr1( t_Start - prev_t_Start_B_addr1 + prev_t_End_B_addr1, i_Start - prev_i_Start_B_addr1 + prev_i_End_B_addr1, j_Start - prev_j_Start_B_addr1 + prev_j_End_B_addr1) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr2 != -1) {
            if ( calAddrB_addr2( t_Start - prev_t_Start_B_addr2 + prev_t_End_B_addr2, i_Start - prev_i_Start_B_addr2 + prev_i_End_B_addr2, j_Start - prev_j_Start_B_addr2 + prev_j_End_B_addr2) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr3 != -1) {
            if ( calAddrB_addr3( t_Start - prev_t_Start_B_addr3 + prev_t_End_B_addr3, i_Start - prev_i_Start_B_addr3 + prev_i_End_B_addr3, j_Start - prev_j_Start_B_addr3 + prev_j_End_B_addr3) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr4 != -1) {
            if ( calAddrB_addr4( t_Start - prev_t_Start_B_addr4 + prev_t_End_B_addr4, i_Start - prev_i_Start_B_addr4 + prev_i_End_B_addr4, j_Start - prev_j_Start_B_addr4 + prev_j_End_B_addr4) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr5 != -1) {
            if ( calAddrB_addr5( t_Start - prev_t_Start_B_addr5 + prev_t_End_B_addr5, i_Start - prev_i_Start_B_addr5 + prev_i_End_B_addr5, j_Start - prev_j_Start_B_addr5 + prev_j_End_B_addr5) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( t, i, j) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB3 = 1;
            if ( t == t_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( t, i, j) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr1 = cnt;
                            prev_t_Start_B_addr1 = t_Start;
                            prev_t_End_B_addr1 = t;
                            prev_i_Start_B_addr1 = i_Start;
                            prev_i_End_B_addr1 = i;
                            prev_j_Start_B_addr1 = j_Start;
                            prev_j_End_B_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( t, i, j) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr2 = cnt;
                            prev_t_Start_B_addr2 = t_Start;
                            prev_t_End_B_addr2 = t;
                            prev_i_Start_B_addr2 = i_Start;
                            prev_i_End_B_addr2 = i;
                            prev_j_Start_B_addr2 = j_Start;
                            prev_j_End_B_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr3( t, i, j) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr3 = cnt;
                            prev_t_Start_B_addr3 = t_Start;
                            prev_t_End_B_addr3 = t;
                            prev_i_Start_B_addr3 = i_Start;
                            prev_i_End_B_addr3 = i;
                            prev_j_Start_B_addr3 = j_Start;
                            prev_j_End_B_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr4( t, i, j) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr4 = cnt;
                            prev_t_Start_B_addr4 = t_Start;
                            prev_t_End_B_addr4 = t;
                            prev_i_Start_B_addr4 = i_Start;
                            prev_i_End_B_addr4 = i;
                            prev_j_Start_B_addr4 = j_Start;
                            prev_j_End_B_addr4 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr5( t, i, j) == calAddrB_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr5 = cnt;
                            prev_t_Start_B_addr5 = t_Start;
                            prev_t_End_B_addr5 = t;
                            prev_i_Start_B_addr5 = i_Start;
                            prev_i_End_B_addr5 = i;
                            prev_j_Start_B_addr5 = j_Start;
                            prev_j_End_B_addr5 = j;
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
void ref_B_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr1 = -1;
    uint64_t prev_t_Start_B_addr1 = -1;
    uint64_t prev_t_End_B_addr1 = -1;
    uint64_t prev_i_Start_B_addr1 = -1;
    uint64_t prev_i_End_B_addr1 = -1;
    uint64_t prev_j_Start_B_addr1 = -1;
    uint64_t prev_j_End_B_addr1 = -1;
    uint64_t prev_cnt_B_addr2 = -1;
    uint64_t prev_t_Start_B_addr2 = -1;
    uint64_t prev_t_End_B_addr2 = -1;
    uint64_t prev_i_Start_B_addr2 = -1;
    uint64_t prev_i_End_B_addr2 = -1;
    uint64_t prev_j_Start_B_addr2 = -1;
    uint64_t prev_j_End_B_addr2 = -1;
    uint64_t prev_cnt_B_addr3 = -1;
    uint64_t prev_t_Start_B_addr3 = -1;
    uint64_t prev_t_End_B_addr3 = -1;
    uint64_t prev_i_Start_B_addr3 = -1;
    uint64_t prev_i_End_B_addr3 = -1;
    uint64_t prev_j_Start_B_addr3 = -1;
    uint64_t prev_j_End_B_addr3 = -1;
    uint64_t prev_cnt_B_addr4 = -1;
    uint64_t prev_t_Start_B_addr4 = -1;
    uint64_t prev_t_End_B_addr4 = -1;
    uint64_t prev_i_Start_B_addr4 = -1;
    uint64_t prev_i_End_B_addr4 = -1;
    uint64_t prev_j_Start_B_addr4 = -1;
    uint64_t prev_j_End_B_addr4 = -1;
    uint64_t prev_cnt_B_addr5 = -1;
    uint64_t prev_t_Start_B_addr5 = -1;
    uint64_t prev_t_End_B_addr5 = -1;
    uint64_t prev_i_Start_B_addr5 = -1;
    uint64_t prev_i_End_B_addr5 = -1;
    uint64_t prev_j_Start_B_addr5 = -1;
    uint64_t prev_j_End_B_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr1 != -1) {
            if ( calAddrB_addr1( t_Start - prev_t_Start_B_addr1 + prev_t_End_B_addr1, i_Start - prev_i_Start_B_addr1 + prev_i_End_B_addr1, j_Start - prev_j_Start_B_addr1 + prev_j_End_B_addr1) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr2 != -1) {
            if ( calAddrB_addr2( t_Start - prev_t_Start_B_addr2 + prev_t_End_B_addr2, i_Start - prev_i_Start_B_addr2 + prev_i_End_B_addr2, j_Start - prev_j_Start_B_addr2 + prev_j_End_B_addr2) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr3 != -1) {
            if ( calAddrB_addr3( t_Start - prev_t_Start_B_addr3 + prev_t_End_B_addr3, i_Start - prev_i_Start_B_addr3 + prev_i_End_B_addr3, j_Start - prev_j_Start_B_addr3 + prev_j_End_B_addr3) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr4 != -1) {
            if ( calAddrB_addr4( t_Start - prev_t_Start_B_addr4 + prev_t_End_B_addr4, i_Start - prev_i_Start_B_addr4 + prev_i_End_B_addr4, j_Start - prev_j_Start_B_addr4 + prev_j_End_B_addr4) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr5 != -1) {
            if ( calAddrB_addr5( t_Start - prev_t_Start_B_addr5 + prev_t_End_B_addr5, i_Start - prev_i_Start_B_addr5 + prev_i_End_B_addr5, j_Start - prev_j_Start_B_addr5 + prev_j_End_B_addr5) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( t, i, j) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB3 = 1;
            if ( t == t_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1023; i++) {
                {
                int jLB4 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( t, i, j) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr1 = cnt;
                            prev_t_Start_B_addr1 = t_Start;
                            prev_t_End_B_addr1 = t;
                            prev_i_Start_B_addr1 = i_Start;
                            prev_i_End_B_addr1 = i;
                            prev_j_Start_B_addr1 = j_Start;
                            prev_j_End_B_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( t, i, j) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr2 = cnt;
                            prev_t_Start_B_addr2 = t_Start;
                            prev_t_End_B_addr2 = t;
                            prev_i_Start_B_addr2 = i_Start;
                            prev_i_End_B_addr2 = i;
                            prev_j_Start_B_addr2 = j_Start;
                            prev_j_End_B_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr3( t, i, j) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr3 = cnt;
                            prev_t_Start_B_addr3 = t_Start;
                            prev_t_End_B_addr3 = t;
                            prev_i_Start_B_addr3 = i_Start;
                            prev_i_End_B_addr3 = i;
                            prev_j_Start_B_addr3 = j_Start;
                            prev_j_End_B_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr4( t, i, j) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr4 = cnt;
                            prev_t_Start_B_addr4 = t_Start;
                            prev_t_End_B_addr4 = t;
                            prev_i_Start_B_addr4 = i_Start;
                            prev_i_End_B_addr4 = i;
                            prev_j_Start_B_addr4 = j_Start;
                            prev_j_End_B_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr5( t, i, j) == calAddrB_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr5 = cnt;
                            prev_t_Start_B_addr5 = t_Start;
                            prev_t_End_B_addr5 = t;
                            prev_i_Start_B_addr5 = i_Start;
                            prev_i_End_B_addr5 = i;
                            prev_j_Start_B_addr5 = j_Start;
                            prev_j_End_B_addr5 = j;
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
    ref_A_addr0();
    ref_A_addr1();
    ref_A_addr2();
    ref_A_addr3();
    ref_A_addr4();
    ref_A_addr5();
    ref_B_addr0();
    ref_B_addr1();
    ref_B_addr2();
    ref_B_addr3();
    ref_B_addr4();
    ref_B_addr5();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: jacobi_2d */ 
