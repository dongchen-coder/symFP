
 /* Start to analysis array index
Array index info
B.addr ((k * 1024) + j)
A.addr ((k * 1024) + i)
B.addr ((i * 1024) + j)
B.addr ((i * 1024) + j)
B.addr ((i * 1024) + j)
B.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %B
double %alpha

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
------k
------Loop Bound: ((i + 1), 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((k * 1024) + i)
--------array access B.addr ((k * 1024) + j)
--------array access B.addr ((i * 1024) + j)
--------array access B.addr ((i * 1024) + j)
------array access B.addr ((i * 1024) + j)
------array access B.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 0 0 1 
Dump stride: 1 1 1 
Dump tree:
----Sample number: 10
------Sample number: 104
--------Sample number: 536
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
int calAddrA_addr0( int i, int j, int k) {
    int result = (((k * 1024) + i)) * 8 / 64;
    return result;
}
int calAddrB_addr0( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrB_addr1( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrB_addr2( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrB_addr3( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrB_addr4( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_k_Start_A_addr0 = -1;
    uint64_t prev_k_End_A_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0, k_Start - prev_k_Start_A_addr0 + prev_k_End_A_addr0) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                {
                int kLB2 = (i + 1);
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            prev_k_Start_A_addr0 = k_Start;
                            prev_k_End_A_addr0 = k;
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
                if (cntStart == true) cnt++;
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
    uint64_t prev_i_Start_B_addr0 = -1;
    uint64_t prev_i_End_B_addr0 = -1;
    uint64_t prev_j_Start_B_addr0 = -1;
    uint64_t prev_j_End_B_addr0 = -1;
    uint64_t prev_k_Start_B_addr0 = -1;
    uint64_t prev_k_End_B_addr0 = -1;
    uint64_t prev_cnt_B_addr1 = -1;
    uint64_t prev_i_Start_B_addr1 = -1;
    uint64_t prev_i_End_B_addr1 = -1;
    uint64_t prev_j_Start_B_addr1 = -1;
    uint64_t prev_j_End_B_addr1 = -1;
    uint64_t prev_k_Start_B_addr1 = -1;
    uint64_t prev_k_End_B_addr1 = -1;
    uint64_t prev_cnt_B_addr2 = -1;
    uint64_t prev_i_Start_B_addr2 = -1;
    uint64_t prev_i_End_B_addr2 = -1;
    uint64_t prev_j_Start_B_addr2 = -1;
    uint64_t prev_j_End_B_addr2 = -1;
    uint64_t prev_k_Start_B_addr2 = -1;
    uint64_t prev_k_End_B_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr0 != -1) {
            if ( calAddrB_addr0( i_Start - prev_i_Start_B_addr0 + prev_i_End_B_addr0, j_Start - prev_j_Start_B_addr0 + prev_j_End_B_addr0, k_Start - prev_k_Start_B_addr0 + prev_k_End_B_addr0) == calAddrB_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr1 != -1) {
            if ( calAddrB_addr1( i_Start - prev_i_Start_B_addr1 + prev_i_End_B_addr1, j_Start - prev_j_Start_B_addr1 + prev_j_End_B_addr1, k_Start - prev_k_Start_B_addr1 + prev_k_End_B_addr1) == calAddrB_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr2 != -1) {
            if ( calAddrB_addr2( i_Start - prev_i_Start_B_addr2 + prev_i_End_B_addr2, j_Start - prev_j_Start_B_addr2 + prev_j_End_B_addr2, k_Start - prev_k_Start_B_addr2 + prev_k_End_B_addr2) == calAddrB_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                {
                int kLB2 = (i + 1);
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( i, j, k) == calAddrB_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr0 = cnt;
                            prev_i_Start_B_addr0 = i_Start;
                            prev_i_End_B_addr0 = i;
                            prev_j_Start_B_addr0 = j_Start;
                            prev_j_End_B_addr0 = j;
                            prev_k_Start_B_addr0 = k_Start;
                            prev_k_End_B_addr0 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( i, j, k) == calAddrB_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr1 = cnt;
                            prev_i_Start_B_addr1 = i_Start;
                            prev_i_End_B_addr1 = i;
                            prev_j_Start_B_addr1 = j_Start;
                            prev_j_End_B_addr1 = j;
                            prev_k_Start_B_addr1 = k_Start;
                            prev_k_End_B_addr1 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( i, j, k) == calAddrB_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr2 = cnt;
                            prev_i_Start_B_addr2 = i_Start;
                            prev_i_End_B_addr2 = i;
                            prev_j_Start_B_addr2 = j_Start;
                            prev_j_End_B_addr2 = j;
                            prev_k_Start_B_addr2 = k_Start;
                            prev_k_End_B_addr2 = k;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr3( i, j) == calAddrB_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr4( i, j) == calAddrB_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
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
void ref_B_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr0 = -1;
    uint64_t prev_i_Start_B_addr0 = -1;
    uint64_t prev_i_End_B_addr0 = -1;
    uint64_t prev_j_Start_B_addr0 = -1;
    uint64_t prev_j_End_B_addr0 = -1;
    uint64_t prev_k_Start_B_addr0 = -1;
    uint64_t prev_k_End_B_addr0 = -1;
    uint64_t prev_cnt_B_addr1 = -1;
    uint64_t prev_i_Start_B_addr1 = -1;
    uint64_t prev_i_End_B_addr1 = -1;
    uint64_t prev_j_Start_B_addr1 = -1;
    uint64_t prev_j_End_B_addr1 = -1;
    uint64_t prev_k_Start_B_addr1 = -1;
    uint64_t prev_k_End_B_addr1 = -1;
    uint64_t prev_cnt_B_addr2 = -1;
    uint64_t prev_i_Start_B_addr2 = -1;
    uint64_t prev_i_End_B_addr2 = -1;
    uint64_t prev_j_Start_B_addr2 = -1;
    uint64_t prev_j_End_B_addr2 = -1;
    uint64_t prev_k_Start_B_addr2 = -1;
    uint64_t prev_k_End_B_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr0 != -1) {
            if ( calAddrB_addr0( i_Start - prev_i_Start_B_addr0 + prev_i_End_B_addr0, j_Start - prev_j_Start_B_addr0 + prev_j_End_B_addr0, k_Start - prev_k_Start_B_addr0 + prev_k_End_B_addr0) == calAddrB_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr1 != -1) {
            if ( calAddrB_addr1( i_Start - prev_i_Start_B_addr1 + prev_i_End_B_addr1, j_Start - prev_j_Start_B_addr1 + prev_j_End_B_addr1, k_Start - prev_k_Start_B_addr1 + prev_k_End_B_addr1) == calAddrB_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr2 != -1) {
            if ( calAddrB_addr2( i_Start - prev_i_Start_B_addr2 + prev_i_End_B_addr2, j_Start - prev_j_Start_B_addr2 + prev_j_End_B_addr2, k_Start - prev_k_Start_B_addr2 + prev_k_End_B_addr2) == calAddrB_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                {
                int kLB2 = (i + 1);
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( i, j, k) == calAddrB_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr0 = cnt;
                            prev_i_Start_B_addr0 = i_Start;
                            prev_i_End_B_addr0 = i;
                            prev_j_Start_B_addr0 = j_Start;
                            prev_j_End_B_addr0 = j;
                            prev_k_Start_B_addr0 = k_Start;
                            prev_k_End_B_addr0 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( i, j, k) == calAddrB_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr1 = cnt;
                            prev_i_Start_B_addr1 = i_Start;
                            prev_i_End_B_addr1 = i;
                            prev_j_Start_B_addr1 = j_Start;
                            prev_j_End_B_addr1 = j;
                            prev_k_Start_B_addr1 = k_Start;
                            prev_k_End_B_addr1 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( i, j, k) == calAddrB_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr2 = cnt;
                            prev_i_Start_B_addr2 = i_Start;
                            prev_i_End_B_addr2 = i;
                            prev_j_Start_B_addr2 = j_Start;
                            prev_j_End_B_addr2 = j;
                            prev_k_Start_B_addr2 = k_Start;
                            prev_k_End_B_addr2 = k;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr3( i, j) == calAddrB_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr4( i, j) == calAddrB_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
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
void ref_B_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr0 = -1;
    uint64_t prev_i_Start_B_addr0 = -1;
    uint64_t prev_i_End_B_addr0 = -1;
    uint64_t prev_j_Start_B_addr0 = -1;
    uint64_t prev_j_End_B_addr0 = -1;
    uint64_t prev_k_Start_B_addr0 = -1;
    uint64_t prev_k_End_B_addr0 = -1;
    uint64_t prev_cnt_B_addr1 = -1;
    uint64_t prev_i_Start_B_addr1 = -1;
    uint64_t prev_i_End_B_addr1 = -1;
    uint64_t prev_j_Start_B_addr1 = -1;
    uint64_t prev_j_End_B_addr1 = -1;
    uint64_t prev_k_Start_B_addr1 = -1;
    uint64_t prev_k_End_B_addr1 = -1;
    uint64_t prev_cnt_B_addr2 = -1;
    uint64_t prev_i_Start_B_addr2 = -1;
    uint64_t prev_i_End_B_addr2 = -1;
    uint64_t prev_j_Start_B_addr2 = -1;
    uint64_t prev_j_End_B_addr2 = -1;
    uint64_t prev_k_Start_B_addr2 = -1;
    uint64_t prev_k_End_B_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr0 != -1) {
            if ( calAddrB_addr0( i_Start - prev_i_Start_B_addr0 + prev_i_End_B_addr0, j_Start - prev_j_Start_B_addr0 + prev_j_End_B_addr0, k_Start - prev_k_Start_B_addr0 + prev_k_End_B_addr0) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr1 != -1) {
            if ( calAddrB_addr1( i_Start - prev_i_Start_B_addr1 + prev_i_End_B_addr1, j_Start - prev_j_Start_B_addr1 + prev_j_End_B_addr1, k_Start - prev_k_Start_B_addr1 + prev_k_End_B_addr1) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr2 != -1) {
            if ( calAddrB_addr2( i_Start - prev_i_Start_B_addr2 + prev_i_End_B_addr2, j_Start - prev_j_Start_B_addr2 + prev_j_End_B_addr2, k_Start - prev_k_Start_B_addr2 + prev_k_End_B_addr2) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                {
                int kLB2 = (i + 1);
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( i, j, k) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr0 = cnt;
                            prev_i_Start_B_addr0 = i_Start;
                            prev_i_End_B_addr0 = i;
                            prev_j_Start_B_addr0 = j_Start;
                            prev_j_End_B_addr0 = j;
                            prev_k_Start_B_addr0 = k_Start;
                            prev_k_End_B_addr0 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( i, j, k) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr1 = cnt;
                            prev_i_Start_B_addr1 = i_Start;
                            prev_i_End_B_addr1 = i;
                            prev_j_Start_B_addr1 = j_Start;
                            prev_j_End_B_addr1 = j;
                            prev_k_Start_B_addr1 = k_Start;
                            prev_k_End_B_addr1 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( i, j, k) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_B_addr2 = cnt;
                            prev_i_Start_B_addr2 = i_Start;
                            prev_i_End_B_addr2 = i;
                            prev_j_Start_B_addr2 = j_Start;
                            prev_j_End_B_addr2 = j;
                            prev_k_Start_B_addr2 = k_Start;
                            prev_k_End_B_addr2 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr3( i, j) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr4( i, j) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
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
void ref_B_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr3 = -1;
    uint64_t prev_i_Start_B_addr3 = -1;
    uint64_t prev_i_End_B_addr3 = -1;
    uint64_t prev_j_Start_B_addr3 = -1;
    uint64_t prev_j_End_B_addr3 = -1;
    uint64_t prev_cnt_B_addr4 = -1;
    uint64_t prev_i_Start_B_addr4 = -1;
    uint64_t prev_i_End_B_addr4 = -1;
    uint64_t prev_j_Start_B_addr4 = -1;
    uint64_t prev_j_End_B_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr3 != -1) {
            if ( calAddrB_addr3( i_Start - prev_i_Start_B_addr3 + prev_i_End_B_addr3, j_Start - prev_j_Start_B_addr3 + prev_j_End_B_addr3) == calAddrB_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr4 != -1) {
            if ( calAddrB_addr4( i_Start - prev_i_Start_B_addr4 + prev_i_End_B_addr4, j_Start - prev_j_Start_B_addr4 + prev_j_End_B_addr4) == calAddrB_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( i, j, k) == calAddrB_addr3(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( i, j, k) == calAddrB_addr3(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( i, j, k) == calAddrB_addr3(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr3( i, j) == calAddrB_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_B_addr3 = cnt;
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
                    if ( calAddrB_addr4( i, j) == calAddrB_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_B_addr4 = cnt;
                        prev_i_Start_B_addr4 = i_Start;
                        prev_i_End_B_addr4 = i;
                        prev_j_Start_B_addr4 = j_Start;
                        prev_j_End_B_addr4 = j;
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
void ref_B_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr3 = -1;
    uint64_t prev_i_Start_B_addr3 = -1;
    uint64_t prev_i_End_B_addr3 = -1;
    uint64_t prev_j_Start_B_addr3 = -1;
    uint64_t prev_j_End_B_addr3 = -1;
    uint64_t prev_cnt_B_addr4 = -1;
    uint64_t prev_i_Start_B_addr4 = -1;
    uint64_t prev_i_End_B_addr4 = -1;
    uint64_t prev_j_Start_B_addr4 = -1;
    uint64_t prev_j_End_B_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr3 != -1) {
            if ( calAddrB_addr3( i_Start - prev_i_Start_B_addr3 + prev_i_End_B_addr3, j_Start - prev_j_Start_B_addr3 + prev_j_End_B_addr3) == calAddrB_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_B_addr4 != -1) {
            if ( calAddrB_addr4( i_Start - prev_i_Start_B_addr4 + prev_i_End_B_addr4, j_Start - prev_j_Start_B_addr4 + prev_j_End_B_addr4) == calAddrB_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_B_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr0( i, j, k) == calAddrB_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr1( i, j, k) == calAddrB_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( i, j, k) == calAddrB_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr3( i, j) == calAddrB_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_B_addr3 = cnt;
                        prev_i_Start_B_addr3 = i_Start;
                        prev_i_End_B_addr3 = i;
                        prev_j_Start_B_addr3 = j_Start;
                        prev_j_End_B_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr4( i, j) == calAddrB_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_B_addr4 = cnt;
                        prev_i_Start_B_addr4 = i_Start;
                        prev_i_End_B_addr4 = i;
                        prev_j_Start_B_addr4 = j_Start;
                        prev_j_End_B_addr4 = j;
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
int main() {
    ref_A_addr0();
    ref_B_addr0();
    ref_B_addr1();
    ref_B_addr2();
    ref_B_addr3();
    ref_B_addr4();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
trmm */ 