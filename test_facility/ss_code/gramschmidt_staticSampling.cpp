
 /* Start to analysis array index
Array index info
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + k)
R.addr ((k * 1024) + k)
A.addr ((i * 1024) + k)
R.addr ((k * 1024) + k)
Q.addr ((i * 1024) + k)
R.addr ((k * 1024) + j)
Q.addr ((i * 1024) + k)
A.addr ((i * 1024) + j)
R.addr ((k * 1024) + j)
R.addr ((k * 1024) + j)
A.addr ((i * 1024) + j)
Q.addr ((i * 1024) + k)
R.addr ((k * 1024) + j)
A.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %R
double* %Q

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--k
--Loop Bound: (0, 1024)
--Loop inc: (k + 1)
--Loop predicate: <
----i
----Loop Bound: (0, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------array access A.addr ((i * 1024) + k)
------array access A.addr ((i * 1024) + k)
----array access R.addr ((k * 1024) + k)
----i
----Loop Bound: (0, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------array access A.addr ((i * 1024) + k)
------array access R.addr ((k * 1024) + k)
------array access Q.addr ((i * 1024) + k)
----j
----Loop Bound: ((k + 1), 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access R.addr ((k * 1024) + j)
------i
------Loop Bound: (0, 1024)
------Loop inc: (i + 1)
------Loop predicate: <
--------array access Q.addr ((i * 1024) + k)
--------array access A.addr ((i * 1024) + j)
--------array access R.addr ((k * 1024) + j)
--------array access R.addr ((k * 1024) + j)
------i
------Loop Bound: (0, 1024)
------Loop inc: (i + 1)
------Loop predicate: <
--------array access A.addr ((i * 1024) + j)
--------array access Q.addr ((i * 1024) + k)
--------array access R.addr ((k * 1024) + j)
--------array access A.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 0 1 
Dump stride: 1 1 
init counter: 0 1 0 
Dump stride: 1 1 1 
init counter: 0 1 0 
Dump stride: 1 1 1 
Dump tree:
----Sample number: 10
------Sample number: 104
------Sample number: 104
------Sample number: 52
--------Sample number: 536
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
int calAddrA_addr0( int k, int i) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrA_addr1( int k, int i) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrR_addr0( int k) {
    int result = (((k * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrA_addr2( int k, int i) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrR_addr1( int k, int i) {
    int result = (((k * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrQ_addr0( int k, int i) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrR_addr2( int k, int j) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrQ_addr1( int k, int j, int i) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrA_addr3( int k, int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrR_addr3( int k, int j, int i) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrR_addr4( int k, int j, int i) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr4( int k, int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrQ_addr2( int k, int j, int i) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrR_addr5( int k, int j, int i) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr5( int k, int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_k_Start_A_addr0 = -1;
    uint64_t prev_k_End_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_k_Start_A_addr1 = -1;
    uint64_t prev_k_End_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( k_Start - prev_k_Start_A_addr0 + prev_k_End_A_addr0, i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0) == calAddrA_addr0(k_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( k_Start - prev_k_Start_A_addr1 + prev_k_End_A_addr1, i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1) == calAddrA_addr0(k_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( k, i) == calAddrA_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr0 = cnt;
                        prev_k_Start_A_addr0 = k_Start;
                        prev_k_End_A_addr0 = k;
                        prev_i_Start_A_addr0 = i_Start;
                        prev_i_End_A_addr0 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( k, i) == calAddrA_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr1 = cnt;
                        prev_k_Start_A_addr1 = k_Start;
                        prev_k_End_A_addr1 = k;
                        prev_i_Start_A_addr1 = i_Start;
                        prev_i_End_A_addr1 = i;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( k, i) == calAddrA_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( k, j, i) == calAddrA_addr0(k_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( k, j, i) == calAddrA_addr0(k_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( k, j, i) == calAddrA_addr0(k_Start, i_Start)) {
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
    uint64_t prev_k_Start_A_addr0 = -1;
    uint64_t prev_k_End_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_k_Start_A_addr1 = -1;
    uint64_t prev_k_End_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( k_Start - prev_k_Start_A_addr0 + prev_k_End_A_addr0, i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0) == calAddrA_addr1(k_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( k_Start - prev_k_Start_A_addr1 + prev_k_End_A_addr1, i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1) == calAddrA_addr1(k_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( k, i) == calAddrA_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr0 = cnt;
                        prev_k_Start_A_addr0 = k_Start;
                        prev_k_End_A_addr0 = k;
                        prev_i_Start_A_addr0 = i_Start;
                        prev_i_End_A_addr0 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( k, i) == calAddrA_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr1 = cnt;
                        prev_k_Start_A_addr1 = k_Start;
                        prev_k_End_A_addr1 = k;
                        prev_i_Start_A_addr1 = i_Start;
                        prev_i_End_A_addr1 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( k, i) == calAddrA_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( k, j, i) == calAddrA_addr1(k_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( k, j, i) == calAddrA_addr1(k_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( k, j, i) == calAddrA_addr1(k_Start, i_Start)) {
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
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_k_Start_A_addr2 = -1;
    uint64_t prev_k_End_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( k_Start - prev_k_Start_A_addr2 + prev_k_End_A_addr2, i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2) == calAddrA_addr2(k_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( k, i) == calAddrA_addr2(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( k, i) == calAddrA_addr2(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            if ( k == k_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( k, i) == calAddrA_addr2(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr2 = cnt;
                        prev_k_Start_A_addr2 = k_Start;
                        prev_k_End_A_addr2 = k;
                        prev_i_Start_A_addr2 = i_Start;
                        prev_i_End_A_addr2 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( k, j, i) == calAddrA_addr2(k_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( k, j, i) == calAddrA_addr2(k_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( k, j, i) == calAddrA_addr2(k_Start, i_Start)) {
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
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_k_Start_A_addr3 = -1;
    uint64_t prev_k_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (k_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (k_Start + 1)) + (k_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( k_Start - prev_k_Start_A_addr3 + prev_k_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3, i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3) == calAddrA_addr3(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( k, i) == calAddrA_addr3(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( k, i) == calAddrA_addr3(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( k, i) == calAddrA_addr3(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            if ( k == k_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int iLB4 = 0;
                if ( k == k_Start && j == j_Start ) {
                    iLB4 = i_Start;
                }
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( k, j, i) == calAddrA_addr3(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_k_Start_A_addr3 = k_Start;
                            prev_k_End_A_addr3 = k;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( k, j, i) == calAddrA_addr3(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( k, j, i) == calAddrA_addr3(k_Start, j_Start, i_Start)) {
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
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_k_Start_A_addr4 = -1;
    uint64_t prev_k_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_cnt_A_addr5 = -1;
    uint64_t prev_k_Start_A_addr5 = -1;
    uint64_t prev_k_End_A_addr5 = -1;
    uint64_t prev_j_Start_A_addr5 = -1;
    uint64_t prev_j_End_A_addr5 = -1;
    uint64_t prev_i_Start_A_addr5 = -1;
    uint64_t prev_i_End_A_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (k_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (k_Start + 1)) + (k_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( k_Start - prev_k_Start_A_addr4 + prev_k_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4, i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4) == calAddrA_addr4(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr5 != -1) {
            if ( calAddrA_addr5( k_Start - prev_k_Start_A_addr5 + prev_k_End_A_addr5, j_Start - prev_j_Start_A_addr5 + prev_j_End_A_addr5, i_Start - prev_i_Start_A_addr5 + prev_i_End_A_addr5) == calAddrA_addr4(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( k, i) == calAddrA_addr4(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( k, i) == calAddrA_addr4(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( k, i) == calAddrA_addr4(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            if ( k == k_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( k, j, i) == calAddrA_addr4(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                {
                int iLB5 = 0;
                if ( k == k_Start && j == j_Start ) {
                    iLB5 = i_Start;
                }
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( k, j, i) == calAddrA_addr4(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr4 = cnt;
                            prev_k_Start_A_addr4 = k_Start;
                            prev_k_End_A_addr4 = k;
                            prev_j_Start_A_addr4 = j_Start;
                            prev_j_End_A_addr4 = j;
                            prev_i_Start_A_addr4 = i_Start;
                            prev_i_End_A_addr4 = i;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( k, j, i) == calAddrA_addr4(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr5 = cnt;
                            prev_k_Start_A_addr5 = k_Start;
                            prev_k_End_A_addr5 = k;
                            prev_j_Start_A_addr5 = j_Start;
                            prev_j_End_A_addr5 = j;
                            prev_i_Start_A_addr5 = i_Start;
                            prev_i_End_A_addr5 = i;
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
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_k_Start_A_addr4 = -1;
    uint64_t prev_k_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_cnt_A_addr5 = -1;
    uint64_t prev_k_Start_A_addr5 = -1;
    uint64_t prev_k_End_A_addr5 = -1;
    uint64_t prev_j_Start_A_addr5 = -1;
    uint64_t prev_j_End_A_addr5 = -1;
    uint64_t prev_i_Start_A_addr5 = -1;
    uint64_t prev_i_End_A_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (k_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (k_Start + 1)) + (k_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( k_Start - prev_k_Start_A_addr4 + prev_k_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4, i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4) == calAddrA_addr5(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr5 != -1) {
            if ( calAddrA_addr5( k_Start - prev_k_Start_A_addr5 + prev_k_End_A_addr5, j_Start - prev_j_Start_A_addr5 + prev_j_End_A_addr5, i_Start - prev_i_Start_A_addr5 + prev_i_End_A_addr5) == calAddrA_addr5(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_A_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( k, i) == calAddrA_addr5(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( k, i) == calAddrA_addr5(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( k, i) == calAddrA_addr5(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            if ( k == k_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( k, j, i) == calAddrA_addr5(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                {
                int iLB5 = 0;
                if ( k == k_Start && j == j_Start ) {
                    iLB5 = i_Start;
                }
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( k, j, i) == calAddrA_addr5(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr4 = cnt;
                            prev_k_Start_A_addr4 = k_Start;
                            prev_k_End_A_addr4 = k;
                            prev_j_Start_A_addr4 = j_Start;
                            prev_j_End_A_addr4 = j;
                            prev_i_Start_A_addr4 = i_Start;
                            prev_i_End_A_addr4 = i;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( k, j, i) == calAddrA_addr5(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr5 = cnt;
                            prev_k_Start_A_addr5 = k_Start;
                            prev_k_End_A_addr5 = k;
                            prev_j_Start_A_addr5 = j_Start;
                            prev_j_End_A_addr5 = j;
                            prev_i_Start_A_addr5 = i_Start;
                            prev_i_End_A_addr5 = i;
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
void ref_Q_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_Q_addr0 = -1;
    uint64_t prev_k_Start_Q_addr0 = -1;
    uint64_t prev_k_End_Q_addr0 = -1;
    uint64_t prev_i_Start_Q_addr0 = -1;
    uint64_t prev_i_End_Q_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_Q_addr0 != -1) {
            if ( calAddrQ_addr0( k_Start - prev_k_Start_Q_addr0 + prev_k_End_Q_addr0, i_Start - prev_i_Start_Q_addr0 + prev_i_End_Q_addr0) == calAddrQ_addr0(k_Start, i_Start)) {
                rtHistoCal(prev_cnt_Q_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            if ( k == k_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrQ_addr0( k, i) == calAddrQ_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_Q_addr0 = cnt;
                        prev_k_Start_Q_addr0 = k_Start;
                        prev_k_End_Q_addr0 = k;
                        prev_i_Start_Q_addr0 = i_Start;
                        prev_i_End_Q_addr0 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            {
            int jLB3 = (k + 1);
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrQ_addr1( k, j, i) == calAddrQ_addr0(k_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrQ_addr2( k, j, i) == calAddrQ_addr0(k_Start, i_Start)) {
                            rtHistoCal(cnt);
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
void ref_Q_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_Q_addr1 = -1;
    uint64_t prev_k_Start_Q_addr1 = -1;
    uint64_t prev_k_End_Q_addr1 = -1;
    uint64_t prev_j_Start_Q_addr1 = -1;
    uint64_t prev_j_End_Q_addr1 = -1;
    uint64_t prev_i_Start_Q_addr1 = -1;
    uint64_t prev_i_End_Q_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (k_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (k_Start + 1)) + (k_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_Q_addr1 != -1) {
            if ( calAddrQ_addr1( k_Start - prev_k_Start_Q_addr1 + prev_k_End_Q_addr1, j_Start - prev_j_Start_Q_addr1 + prev_j_End_Q_addr1, i_Start - prev_i_Start_Q_addr1 + prev_i_End_Q_addr1) == calAddrQ_addr1(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_Q_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrQ_addr0( k, i) == calAddrQ_addr1(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int jLB3 = (k + 1);
            if ( k == k_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int iLB4 = 0;
                if ( k == k_Start && j == j_Start ) {
                    iLB4 = i_Start;
                }
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrQ_addr1( k, j, i) == calAddrQ_addr1(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_Q_addr1 = cnt;
                            prev_k_Start_Q_addr1 = k_Start;
                            prev_k_End_Q_addr1 = k;
                            prev_j_Start_Q_addr1 = j_Start;
                            prev_j_End_Q_addr1 = j;
                            prev_i_Start_Q_addr1 = i_Start;
                            prev_i_End_Q_addr1 = i;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrQ_addr2( k, j, i) == calAddrQ_addr1(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
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
void ref_Q_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_Q_addr2 = -1;
    uint64_t prev_k_Start_Q_addr2 = -1;
    uint64_t prev_k_End_Q_addr2 = -1;
    uint64_t prev_j_Start_Q_addr2 = -1;
    uint64_t prev_j_End_Q_addr2 = -1;
    uint64_t prev_i_Start_Q_addr2 = -1;
    uint64_t prev_i_End_Q_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (k_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (k_Start + 1)) + (k_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_Q_addr2 != -1) {
            if ( calAddrQ_addr2( k_Start - prev_k_Start_Q_addr2 + prev_k_End_Q_addr2, j_Start - prev_j_Start_Q_addr2 + prev_j_End_Q_addr2, i_Start - prev_i_Start_Q_addr2 + prev_i_End_Q_addr2) == calAddrQ_addr2(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_Q_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrQ_addr0( k, i) == calAddrQ_addr2(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int jLB3 = (k + 1);
            if ( k == k_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrQ_addr1( k, j, i) == calAddrQ_addr2(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                {
                int iLB5 = 0;
                if ( k == k_Start && j == j_Start ) {
                    iLB5 = i_Start;
                }
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrQ_addr2( k, j, i) == calAddrQ_addr2(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_Q_addr2 = cnt;
                            prev_k_Start_Q_addr2 = k_Start;
                            prev_k_End_Q_addr2 = k;
                            prev_j_Start_Q_addr2 = j_Start;
                            prev_j_End_Q_addr2 = j;
                            prev_i_Start_Q_addr2 = i_Start;
                            prev_i_End_Q_addr2 = i;
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
void ref_R_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_R_addr0 = -1;
    uint64_t prev_k_Start_R_addr0 = -1;
    uint64_t prev_k_End_R_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_R_addr0 != -1) {
            if ( calAddrR_addr0( k_Start - prev_k_Start_R_addr0 + prev_k_End_R_addr0) == calAddrR_addr0(k_Start)) {
                rtHistoCal(prev_cnt_R_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrR_addr0( k) == calAddrR_addr0(k_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_R_addr0 = cnt;
                    prev_k_Start_R_addr0 = k_Start;
                    prev_k_End_R_addr0 = k;
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr1( k, i) == calAddrR_addr0(k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr2( k, j) == calAddrR_addr0(k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr3( k, j, i) == calAddrR_addr0(k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr4( k, j, i) == calAddrR_addr0(k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr5( k, j, i) == calAddrR_addr0(k_Start)) {
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
void ref_R_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_R_addr1 = -1;
    uint64_t prev_k_Start_R_addr1 = -1;
    uint64_t prev_k_End_R_addr1 = -1;
    uint64_t prev_i_Start_R_addr1 = -1;
    uint64_t prev_i_End_R_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_R_addr1 != -1) {
            if ( calAddrR_addr1( k_Start - prev_k_Start_R_addr1 + prev_k_End_R_addr1, i_Start - prev_i_Start_R_addr1 + prev_i_End_R_addr1) == calAddrR_addr1(k_Start, i_Start)) {
                rtHistoCal(prev_cnt_R_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrR_addr0( k) == calAddrR_addr1(k_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB2 = 0;
            if ( k == k_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr1( k, i) == calAddrR_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_R_addr1 = cnt;
                        prev_k_Start_R_addr1 = k_Start;
                        prev_k_End_R_addr1 = k;
                        prev_i_Start_R_addr1 = i_Start;
                        prev_i_End_R_addr1 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr2( k, j) == calAddrR_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr3( k, j, i) == calAddrR_addr1(k_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr4( k, j, i) == calAddrR_addr1(k_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr5( k, j, i) == calAddrR_addr1(k_Start, i_Start)) {
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
void ref_R_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_R_addr2 = -1;
    uint64_t prev_k_Start_R_addr2 = -1;
    uint64_t prev_k_End_R_addr2 = -1;
    uint64_t prev_j_Start_R_addr2 = -1;
    uint64_t prev_j_End_R_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (k_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (k_Start + 1)) + (k_Start + 1);
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_R_addr2 != -1) {
            if ( calAddrR_addr2( k_Start - prev_k_Start_R_addr2 + prev_k_End_R_addr2, j_Start - prev_j_Start_R_addr2 + prev_j_End_R_addr2) == calAddrR_addr2(k_Start, j_Start)) {
                rtHistoCal(prev_cnt_R_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrR_addr0( k) == calAddrR_addr2(k_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr1( k, i) == calAddrR_addr2(k_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            if ( k == k_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr2( k, j) == calAddrR_addr2(k_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_R_addr2 = cnt;
                        prev_k_Start_R_addr2 = k_Start;
                        prev_k_End_R_addr2 = k;
                        prev_j_Start_R_addr2 = j_Start;
                        prev_j_End_R_addr2 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr3( k, j, i) == calAddrR_addr2(k_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr4( k, j, i) == calAddrR_addr2(k_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr5( k, j, i) == calAddrR_addr2(k_Start, j_Start)) {
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
void ref_R_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_R_addr3 = -1;
    uint64_t prev_k_Start_R_addr3 = -1;
    uint64_t prev_k_End_R_addr3 = -1;
    uint64_t prev_j_Start_R_addr3 = -1;
    uint64_t prev_j_End_R_addr3 = -1;
    uint64_t prev_i_Start_R_addr3 = -1;
    uint64_t prev_i_End_R_addr3 = -1;
    uint64_t prev_cnt_R_addr4 = -1;
    uint64_t prev_k_Start_R_addr4 = -1;
    uint64_t prev_k_End_R_addr4 = -1;
    uint64_t prev_j_Start_R_addr4 = -1;
    uint64_t prev_j_End_R_addr4 = -1;
    uint64_t prev_i_Start_R_addr4 = -1;
    uint64_t prev_i_End_R_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (k_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (k_Start + 1)) + (k_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_R_addr3 != -1) {
            if ( calAddrR_addr3( k_Start - prev_k_Start_R_addr3 + prev_k_End_R_addr3, j_Start - prev_j_Start_R_addr3 + prev_j_End_R_addr3, i_Start - prev_i_Start_R_addr3 + prev_i_End_R_addr3) == calAddrR_addr3(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_R_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_R_addr4 != -1) {
            if ( calAddrR_addr4( k_Start - prev_k_Start_R_addr4 + prev_k_End_R_addr4, j_Start - prev_j_Start_R_addr4 + prev_j_End_R_addr4, i_Start - prev_i_Start_R_addr4 + prev_i_End_R_addr4) == calAddrR_addr3(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_R_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrR_addr0( k) == calAddrR_addr3(k_Start, j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr1( k, i) == calAddrR_addr3(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            if ( k == k_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr2( k, j) == calAddrR_addr3(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB4 = 0;
                if ( k == k_Start && j == j_Start ) {
                    iLB4 = i_Start;
                }
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr3( k, j, i) == calAddrR_addr3(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_R_addr3 = cnt;
                            prev_k_Start_R_addr3 = k_Start;
                            prev_k_End_R_addr3 = k;
                            prev_j_Start_R_addr3 = j_Start;
                            prev_j_End_R_addr3 = j;
                            prev_i_Start_R_addr3 = i_Start;
                            prev_i_End_R_addr3 = i;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr4( k, j, i) == calAddrR_addr3(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_R_addr4 = cnt;
                            prev_k_Start_R_addr4 = k_Start;
                            prev_k_End_R_addr4 = k;
                            prev_j_Start_R_addr4 = j_Start;
                            prev_j_End_R_addr4 = j;
                            prev_i_Start_R_addr4 = i_Start;
                            prev_i_End_R_addr4 = i;
                            goto EndSample;
                        }
                    }
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr5( k, j, i) == calAddrR_addr3(k_Start, j_Start, i_Start)) {
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
void ref_R_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_R_addr3 = -1;
    uint64_t prev_k_Start_R_addr3 = -1;
    uint64_t prev_k_End_R_addr3 = -1;
    uint64_t prev_j_Start_R_addr3 = -1;
    uint64_t prev_j_End_R_addr3 = -1;
    uint64_t prev_i_Start_R_addr3 = -1;
    uint64_t prev_i_End_R_addr3 = -1;
    uint64_t prev_cnt_R_addr4 = -1;
    uint64_t prev_k_Start_R_addr4 = -1;
    uint64_t prev_k_End_R_addr4 = -1;
    uint64_t prev_j_Start_R_addr4 = -1;
    uint64_t prev_j_End_R_addr4 = -1;
    uint64_t prev_i_Start_R_addr4 = -1;
    uint64_t prev_i_End_R_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (k_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (k_Start + 1)) + (k_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_R_addr3 != -1) {
            if ( calAddrR_addr3( k_Start - prev_k_Start_R_addr3 + prev_k_End_R_addr3, j_Start - prev_j_Start_R_addr3 + prev_j_End_R_addr3, i_Start - prev_i_Start_R_addr3 + prev_i_End_R_addr3) == calAddrR_addr4(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_R_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_R_addr4 != -1) {
            if ( calAddrR_addr4( k_Start - prev_k_Start_R_addr4 + prev_k_End_R_addr4, j_Start - prev_j_Start_R_addr4 + prev_j_End_R_addr4, i_Start - prev_i_Start_R_addr4 + prev_i_End_R_addr4) == calAddrR_addr4(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_R_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrR_addr0( k) == calAddrR_addr4(k_Start, j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr1( k, i) == calAddrR_addr4(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            if ( k == k_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr2( k, j) == calAddrR_addr4(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB4 = 0;
                if ( k == k_Start && j == j_Start ) {
                    iLB4 = i_Start;
                }
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr3( k, j, i) == calAddrR_addr4(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_R_addr3 = cnt;
                            prev_k_Start_R_addr3 = k_Start;
                            prev_k_End_R_addr3 = k;
                            prev_j_Start_R_addr3 = j_Start;
                            prev_j_End_R_addr3 = j;
                            prev_i_Start_R_addr3 = i_Start;
                            prev_i_End_R_addr3 = i;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr4( k, j, i) == calAddrR_addr4(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_R_addr4 = cnt;
                            prev_k_Start_R_addr4 = k_Start;
                            prev_k_End_R_addr4 = k;
                            prev_j_Start_R_addr4 = j_Start;
                            prev_j_End_R_addr4 = j;
                            prev_i_Start_R_addr4 = i_Start;
                            prev_i_End_R_addr4 = i;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                {
                int iLB5 = 0;
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr5( k, j, i) == calAddrR_addr4(k_Start, j_Start, i_Start)) {
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
void ref_R_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_R_addr5 = -1;
    uint64_t prev_k_Start_R_addr5 = -1;
    uint64_t prev_k_End_R_addr5 = -1;
    uint64_t prev_j_Start_R_addr5 = -1;
    uint64_t prev_j_End_R_addr5 = -1;
    uint64_t prev_i_Start_R_addr5 = -1;
    uint64_t prev_i_End_R_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - (k_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (k_Start + 1)) + (k_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_R_addr5 != -1) {
            if ( calAddrR_addr5( k_Start - prev_k_Start_R_addr5 + prev_k_End_R_addr5, j_Start - prev_j_Start_R_addr5 + prev_j_End_R_addr5, i_Start - prev_i_Start_R_addr5 + prev_i_End_R_addr5) == calAddrR_addr5(k_Start, j_Start, i_Start)) {
                rtHistoCal(prev_cnt_R_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrR_addr0( k) == calAddrR_addr5(k_Start, j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr1( k, i) == calAddrR_addr5(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int jLB3 = (k + 1);
            if ( k == k_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrR_addr2( k, j) == calAddrR_addr5(k_Start, j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB4 = 0;
                for ( int i = iLB4; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr3( k, j, i) == calAddrR_addr5(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr4( k, j, i) == calAddrR_addr5(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                {
                int iLB5 = 0;
                if ( k == k_Start && j == j_Start ) {
                    iLB5 = i_Start;
                }
                for ( int i = iLB5; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrR_addr5( k, j, i) == calAddrR_addr5(k_Start, j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_R_addr5 = cnt;
                            prev_k_Start_R_addr5 = k_Start;
                            prev_k_End_R_addr5 = k;
                            prev_j_Start_R_addr5 = j_Start;
                            prev_j_End_R_addr5 = j;
                            prev_i_Start_R_addr5 = i_Start;
                            prev_i_End_R_addr5 = i;
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
    ref_Q_addr0();
    ref_Q_addr1();
    ref_Q_addr2();
    ref_R_addr0();
    ref_R_addr1();
    ref_R_addr2();
    ref_R_addr3();
    ref_R_addr4();
    ref_R_addr5();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
gramschmidt */ 
