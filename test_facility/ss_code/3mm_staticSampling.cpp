
 /* Start to analysis array index
Array index info
E.addr ((i * 256) + j)
A.addr ((i * 256) + k)
B.addr ((k * 256) + j)
E.addr ((i * 256) + j)
E.addr ((i * 256) + j)
F.addr ((i * 256) + j)
C.addr ((i * 256) + k)
D.addr ((k * 256) + j)
F.addr ((i * 256) + j)
F.addr ((i * 256) + j)
G.addr ((i * 256) + j)
E.addr ((i * 256) + k)
F.addr ((k * 256) + j)
G.addr ((i * 256) + j)
G.addr ((i * 256) + j)

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
--Loop Bound: (0, 256)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 256)
----Loop inc: (j + 1)
----Loop predicate: <
------array access E.addr ((i * 256) + j)
------k
------Loop Bound: (0, 256)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 256) + k)
--------array access B.addr ((k * 256) + j)
--------array access E.addr ((i * 256) + j)
--------array access E.addr ((i * 256) + j)
--i
--Loop Bound: (0, 256)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 256)
----Loop inc: (j + 1)
----Loop predicate: <
------array access F.addr ((i * 256) + j)
------k
------Loop Bound: (0, 256)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access C.addr ((i * 256) + k)
--------array access D.addr ((k * 256) + j)
--------array access F.addr ((i * 256) + j)
--------array access F.addr ((i * 256) + j)
--i
--Loop Bound: (0, 256)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 256)
----Loop inc: (j + 1)
----Loop predicate: <
------array access G.addr ((i * 256) + j)
------k
------Loop Bound: (0, 256)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access E.addr ((i * 256) + k)
--------array access F.addr ((k * 256) + j)
--------array access G.addr ((i * 256) + j)
--------array access G.addr ((i * 256) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 2
------Sample number: 6
--------Sample number: 16
----Sample number: 2
------Sample number: 6
--------Sample number: 16
----Sample number: 2
------Sample number: 6
--------Sample number: 16
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
int calAddrE_addr0( int i, int j) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr0( int i, int j, int k) {
    int result = (((i * 256) + k)) * 8 / 64;
    return result;
}
int calAddrB_addr0( int i, int j, int k) {
    int result = (((k * 256) + j)) * 8 / 64;
    return result;
}
int calAddrE_addr1( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
int calAddrE_addr2( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
int calAddrF_addr0( int i, int j) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
int calAddrC_addr0( int i, int j, int k) {
    int result = (((i * 256) + k)) * 8 / 64;
    return result;
}
int calAddrD_addr0( int i, int j, int k) {
    int result = (((k * 256) + j)) * 8 / 64;
    return result;
}
int calAddrF_addr1( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
int calAddrF_addr2( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
int calAddrG_addr0( int i, int j) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
int calAddrE_addr3( int i, int j, int k) {
    int result = (((i * 256) + k)) * 8 / 64;
    return result;
}
int calAddrF_addr3( int i, int j, int k) {
    int result = (((k * 256) + j)) * 8 / 64;
    return result;
}
int calAddrG_addr1( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
int calAddrG_addr2( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
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
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
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
        for ( int i = iLB0; i < 256; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 256; k++) {
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
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 256; k++) {
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
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
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
void ref_B_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr0 = -1;
    uint64_t prev_i_Start_B_addr0 = -1;
    uint64_t prev_i_End_B_addr0 = -1;
    uint64_t prev_j_Start_B_addr0 = -1;
    uint64_t prev_j_End_B_addr0 = -1;
    uint64_t prev_k_Start_B_addr0 = -1;
    uint64_t prev_k_End_B_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
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
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 256; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 256; k++) {
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
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 256; k++) {
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
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
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
void ref_C_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_C_addr0 = -1;
    uint64_t prev_i_Start_C_addr0 = -1;
    uint64_t prev_i_End_C_addr0 = -1;
    uint64_t prev_j_Start_C_addr0 = -1;
    uint64_t prev_j_End_C_addr0 = -1;
    uint64_t prev_k_Start_C_addr0 = -1;
    uint64_t prev_k_End_C_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_C_addr0 != -1) {
            if ( calAddrC_addr0( i_Start - prev_i_Start_C_addr0 + prev_i_End_C_addr0, j_Start - prev_j_Start_C_addr0 + prev_j_End_C_addr0, k_Start - prev_k_Start_C_addr0 + prev_k_End_C_addr0) == calAddrC_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_C_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 256; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrC_addr0( i, j, k) == calAddrC_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_C_addr0 = cnt;
                            prev_i_Start_C_addr0 = i_Start;
                            prev_i_End_C_addr0 = i;
                            prev_j_Start_C_addr0 = j_Start;
                            prev_j_End_C_addr0 = j;
                            prev_k_Start_C_addr0 = k_Start;
                            prev_k_End_C_addr0 = k;
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
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
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
void ref_D_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_D_addr0 = -1;
    uint64_t prev_i_Start_D_addr0 = -1;
    uint64_t prev_i_End_D_addr0 = -1;
    uint64_t prev_j_Start_D_addr0 = -1;
    uint64_t prev_j_End_D_addr0 = -1;
    uint64_t prev_k_Start_D_addr0 = -1;
    uint64_t prev_k_End_D_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_D_addr0 != -1) {
            if ( calAddrD_addr0( i_Start - prev_i_Start_D_addr0 + prev_i_End_D_addr0, j_Start - prev_j_Start_D_addr0 + prev_j_End_D_addr0, k_Start - prev_k_Start_D_addr0 + prev_k_End_D_addr0) == calAddrD_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_D_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr0( i, j, k) == calAddrD_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_D_addr0 = cnt;
                            prev_i_Start_D_addr0 = i_Start;
                            prev_i_End_D_addr0 = i;
                            prev_j_Start_D_addr0 = j_Start;
                            prev_j_End_D_addr0 = j;
                            prev_k_Start_D_addr0 = k_Start;
                            prev_k_End_D_addr0 = k;
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
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
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
void ref_E_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_E_addr0 = -1;
    uint64_t prev_i_Start_E_addr0 = -1;
    uint64_t prev_i_End_E_addr0 = -1;
    uint64_t prev_j_Start_E_addr0 = -1;
    uint64_t prev_j_End_E_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_E_addr0 != -1) {
            if ( calAddrE_addr0( i_Start - prev_i_Start_E_addr0 + prev_i_End_E_addr0, j_Start - prev_j_Start_E_addr0 + prev_j_End_E_addr0) == calAddrE_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_E_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 256; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 256; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrE_addr0( i, j) == calAddrE_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_E_addr0 = cnt;
                        prev_i_Start_E_addr0 = i_Start;
                        prev_i_End_E_addr0 = i;
                        prev_j_Start_E_addr0 = j_Start;
                        prev_j_End_E_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr1( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr2( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
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
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 256; k++) {
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
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
void ref_E_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_E_addr1 = -1;
    uint64_t prev_i_Start_E_addr1 = -1;
    uint64_t prev_i_End_E_addr1 = -1;
    uint64_t prev_j_Start_E_addr1 = -1;
    uint64_t prev_j_End_E_addr1 = -1;
    uint64_t prev_k_Start_E_addr1 = -1;
    uint64_t prev_k_End_E_addr1 = -1;
    uint64_t prev_cnt_E_addr2 = -1;
    uint64_t prev_i_Start_E_addr2 = -1;
    uint64_t prev_i_End_E_addr2 = -1;
    uint64_t prev_j_Start_E_addr2 = -1;
    uint64_t prev_j_End_E_addr2 = -1;
    uint64_t prev_k_Start_E_addr2 = -1;
    uint64_t prev_k_End_E_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_E_addr1 != -1) {
            if ( calAddrE_addr1( i_Start - prev_i_Start_E_addr1 + prev_i_End_E_addr1, j_Start - prev_j_Start_E_addr1 + prev_j_End_E_addr1, k_Start - prev_k_Start_E_addr1 + prev_k_End_E_addr1) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_E_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_E_addr2 != -1) {
            if ( calAddrE_addr2( i_Start - prev_i_Start_E_addr2 + prev_i_End_E_addr2, j_Start - prev_j_Start_E_addr2 + prev_j_End_E_addr2, k_Start - prev_k_Start_E_addr2 + prev_k_End_E_addr2) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_E_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 256; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 256; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrE_addr0( i, j) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr1( i, j, k) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_E_addr1 = cnt;
                            prev_i_Start_E_addr1 = i_Start;
                            prev_i_End_E_addr1 = i;
                            prev_j_Start_E_addr1 = j_Start;
                            prev_j_End_E_addr1 = j;
                            prev_k_Start_E_addr1 = k_Start;
                            prev_k_End_E_addr1 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr2( i, j, k) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_E_addr2 = cnt;
                            prev_i_Start_E_addr2 = i_Start;
                            prev_i_End_E_addr2 = i;
                            prev_j_Start_E_addr2 = j_Start;
                            prev_j_End_E_addr2 = j;
                            prev_k_Start_E_addr2 = k_Start;
                            prev_k_End_E_addr2 = k;
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
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 256; k++) {
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
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
void ref_E_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_E_addr1 = -1;
    uint64_t prev_i_Start_E_addr1 = -1;
    uint64_t prev_i_End_E_addr1 = -1;
    uint64_t prev_j_Start_E_addr1 = -1;
    uint64_t prev_j_End_E_addr1 = -1;
    uint64_t prev_k_Start_E_addr1 = -1;
    uint64_t prev_k_End_E_addr1 = -1;
    uint64_t prev_cnt_E_addr2 = -1;
    uint64_t prev_i_Start_E_addr2 = -1;
    uint64_t prev_i_End_E_addr2 = -1;
    uint64_t prev_j_Start_E_addr2 = -1;
    uint64_t prev_j_End_E_addr2 = -1;
    uint64_t prev_k_Start_E_addr2 = -1;
    uint64_t prev_k_End_E_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_E_addr1 != -1) {
            if ( calAddrE_addr1( i_Start - prev_i_Start_E_addr1 + prev_i_End_E_addr1, j_Start - prev_j_Start_E_addr1 + prev_j_End_E_addr1, k_Start - prev_k_Start_E_addr1 + prev_k_End_E_addr1) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_E_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_E_addr2 != -1) {
            if ( calAddrE_addr2( i_Start - prev_i_Start_E_addr2 + prev_i_End_E_addr2, j_Start - prev_j_Start_E_addr2 + prev_j_End_E_addr2, k_Start - prev_k_Start_E_addr2 + prev_k_End_E_addr2) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_E_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 256; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 256; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrE_addr0( i, j) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr1( i, j, k) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_E_addr1 = cnt;
                            prev_i_Start_E_addr1 = i_Start;
                            prev_i_End_E_addr1 = i;
                            prev_j_Start_E_addr1 = j_Start;
                            prev_j_End_E_addr1 = j;
                            prev_k_Start_E_addr1 = k_Start;
                            prev_k_End_E_addr1 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr2( i, j, k) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_E_addr2 = cnt;
                            prev_i_Start_E_addr2 = i_Start;
                            prev_i_End_E_addr2 = i;
                            prev_j_Start_E_addr2 = j_Start;
                            prev_j_End_E_addr2 = j;
                            prev_k_Start_E_addr2 = k_Start;
                            prev_k_End_E_addr2 = k;
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
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 256; k++) {
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
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
void ref_E_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_E_addr3 = -1;
    uint64_t prev_i_Start_E_addr3 = -1;
    uint64_t prev_i_End_E_addr3 = -1;
    uint64_t prev_j_Start_E_addr3 = -1;
    uint64_t prev_j_End_E_addr3 = -1;
    uint64_t prev_k_Start_E_addr3 = -1;
    uint64_t prev_k_End_E_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_E_addr3 != -1) {
            if ( calAddrE_addr3( i_Start - prev_i_Start_E_addr3 + prev_i_End_E_addr3, j_Start - prev_j_Start_E_addr3 + prev_j_End_E_addr3, k_Start - prev_k_Start_E_addr3 + prev_k_End_E_addr3) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_E_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_E_addr3 = cnt;
                            prev_i_Start_E_addr3 = i_Start;
                            prev_i_End_E_addr3 = i;
                            prev_j_Start_E_addr3 = j_Start;
                            prev_j_End_E_addr3 = j;
                            prev_k_Start_E_addr3 = k_Start;
                            prev_k_End_E_addr3 = k;
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
void ref_F_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_F_addr0 = -1;
    uint64_t prev_i_Start_F_addr0 = -1;
    uint64_t prev_i_End_F_addr0 = -1;
    uint64_t prev_j_Start_F_addr0 = -1;
    uint64_t prev_j_End_F_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_F_addr0 != -1) {
            if ( calAddrF_addr0( i_Start - prev_i_Start_F_addr0 + prev_i_End_F_addr0, j_Start - prev_j_Start_F_addr0 + prev_j_End_F_addr0) == calAddrF_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_F_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrF_addr0( i, j) == calAddrF_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_F_addr0 = cnt;
                        prev_i_Start_F_addr0 = i_Start;
                        prev_i_End_F_addr0 = i;
                        prev_j_Start_F_addr0 = j_Start;
                        prev_j_End_F_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr1( i, j, k) == calAddrF_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr2( i, j, k) == calAddrF_addr0(i_Start, j_Start)) {
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
        {
        int iLB6 = 0;
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr3( i, j, k) == calAddrF_addr0(i_Start, j_Start)) {
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
void ref_F_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_F_addr1 = -1;
    uint64_t prev_i_Start_F_addr1 = -1;
    uint64_t prev_i_End_F_addr1 = -1;
    uint64_t prev_j_Start_F_addr1 = -1;
    uint64_t prev_j_End_F_addr1 = -1;
    uint64_t prev_k_Start_F_addr1 = -1;
    uint64_t prev_k_End_F_addr1 = -1;
    uint64_t prev_cnt_F_addr2 = -1;
    uint64_t prev_i_Start_F_addr2 = -1;
    uint64_t prev_i_End_F_addr2 = -1;
    uint64_t prev_j_Start_F_addr2 = -1;
    uint64_t prev_j_End_F_addr2 = -1;
    uint64_t prev_k_Start_F_addr2 = -1;
    uint64_t prev_k_End_F_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_F_addr1 != -1) {
            if ( calAddrF_addr1( i_Start - prev_i_Start_F_addr1 + prev_i_End_F_addr1, j_Start - prev_j_Start_F_addr1 + prev_j_End_F_addr1, k_Start - prev_k_Start_F_addr1 + prev_k_End_F_addr1) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_F_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_F_addr2 != -1) {
            if ( calAddrF_addr2( i_Start - prev_i_Start_F_addr2 + prev_i_End_F_addr2, j_Start - prev_j_Start_F_addr2 + prev_j_End_F_addr2, k_Start - prev_k_Start_F_addr2 + prev_k_End_F_addr2) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_F_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrF_addr0( i, j) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr1( i, j, k) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_F_addr1 = cnt;
                            prev_i_Start_F_addr1 = i_Start;
                            prev_i_End_F_addr1 = i;
                            prev_j_Start_F_addr1 = j_Start;
                            prev_j_End_F_addr1 = j;
                            prev_k_Start_F_addr1 = k_Start;
                            prev_k_End_F_addr1 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr2( i, j, k) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_F_addr2 = cnt;
                            prev_i_Start_F_addr2 = i_Start;
                            prev_i_End_F_addr2 = i;
                            prev_j_Start_F_addr2 = j_Start;
                            prev_j_End_F_addr2 = j;
                            prev_k_Start_F_addr2 = k_Start;
                            prev_k_End_F_addr2 = k;
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
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr3( i, j, k) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
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
void ref_F_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_F_addr1 = -1;
    uint64_t prev_i_Start_F_addr1 = -1;
    uint64_t prev_i_End_F_addr1 = -1;
    uint64_t prev_j_Start_F_addr1 = -1;
    uint64_t prev_j_End_F_addr1 = -1;
    uint64_t prev_k_Start_F_addr1 = -1;
    uint64_t prev_k_End_F_addr1 = -1;
    uint64_t prev_cnt_F_addr2 = -1;
    uint64_t prev_i_Start_F_addr2 = -1;
    uint64_t prev_i_End_F_addr2 = -1;
    uint64_t prev_j_Start_F_addr2 = -1;
    uint64_t prev_j_End_F_addr2 = -1;
    uint64_t prev_k_Start_F_addr2 = -1;
    uint64_t prev_k_End_F_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_F_addr1 != -1) {
            if ( calAddrF_addr1( i_Start - prev_i_Start_F_addr1 + prev_i_End_F_addr1, j_Start - prev_j_Start_F_addr1 + prev_j_End_F_addr1, k_Start - prev_k_Start_F_addr1 + prev_k_End_F_addr1) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_F_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_F_addr2 != -1) {
            if ( calAddrF_addr2( i_Start - prev_i_Start_F_addr2 + prev_i_End_F_addr2, j_Start - prev_j_Start_F_addr2 + prev_j_End_F_addr2, k_Start - prev_k_Start_F_addr2 + prev_k_End_F_addr2) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_F_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 256; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 256; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrF_addr0( i, j) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr1( i, j, k) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_F_addr1 = cnt;
                            prev_i_Start_F_addr1 = i_Start;
                            prev_i_End_F_addr1 = i;
                            prev_j_Start_F_addr1 = j_Start;
                            prev_j_End_F_addr1 = j;
                            prev_k_Start_F_addr1 = k_Start;
                            prev_k_End_F_addr1 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr2( i, j, k) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_F_addr2 = cnt;
                            prev_i_Start_F_addr2 = i_Start;
                            prev_i_End_F_addr2 = i;
                            prev_j_Start_F_addr2 = j_Start;
                            prev_j_End_F_addr2 = j;
                            prev_k_Start_F_addr2 = k_Start;
                            prev_k_End_F_addr2 = k;
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
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr3( i, j, k) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
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
void ref_F_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_F_addr3 = -1;
    uint64_t prev_i_Start_F_addr3 = -1;
    uint64_t prev_i_End_F_addr3 = -1;
    uint64_t prev_j_Start_F_addr3 = -1;
    uint64_t prev_j_End_F_addr3 = -1;
    uint64_t prev_k_Start_F_addr3 = -1;
    uint64_t prev_k_End_F_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_F_addr3 != -1) {
            if ( calAddrF_addr3( i_Start - prev_i_Start_F_addr3 + prev_i_End_F_addr3, j_Start - prev_j_Start_F_addr3 + prev_j_End_F_addr3, k_Start - prev_k_Start_F_addr3 + prev_k_End_F_addr3) == calAddrF_addr3(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_F_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB8 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr3( i, j, k) == calAddrF_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_F_addr3 = cnt;
                            prev_i_Start_F_addr3 = i_Start;
                            prev_i_End_F_addr3 = i;
                            prev_j_Start_F_addr3 = j_Start;
                            prev_j_End_F_addr3 = j;
                            prev_k_Start_F_addr3 = k_Start;
                            prev_k_End_F_addr3 = k;
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
void ref_G_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_G_addr0 = -1;
    uint64_t prev_i_Start_G_addr0 = -1;
    uint64_t prev_i_End_G_addr0 = -1;
    uint64_t prev_j_Start_G_addr0 = -1;
    uint64_t prev_j_End_G_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 6;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_G_addr0 != -1) {
            if ( calAddrG_addr0( i_Start - prev_i_Start_G_addr0 + prev_i_End_G_addr0, j_Start - prev_j_Start_G_addr0 + prev_j_End_G_addr0) == calAddrG_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_G_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrG_addr0( i, j) == calAddrG_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_G_addr0 = cnt;
                        prev_i_Start_G_addr0 = i_Start;
                        prev_i_End_G_addr0 = i;
                        prev_j_Start_G_addr0 = j_Start;
                        prev_j_End_G_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB8 = 0;
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr1( i, j, k) == calAddrG_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr2( i, j, k) == calAddrG_addr0(i_Start, j_Start)) {
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
void ref_G_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_G_addr1 = -1;
    uint64_t prev_i_Start_G_addr1 = -1;
    uint64_t prev_i_End_G_addr1 = -1;
    uint64_t prev_j_Start_G_addr1 = -1;
    uint64_t prev_j_End_G_addr1 = -1;
    uint64_t prev_k_Start_G_addr1 = -1;
    uint64_t prev_k_End_G_addr1 = -1;
    uint64_t prev_cnt_G_addr2 = -1;
    uint64_t prev_i_Start_G_addr2 = -1;
    uint64_t prev_i_End_G_addr2 = -1;
    uint64_t prev_j_Start_G_addr2 = -1;
    uint64_t prev_j_End_G_addr2 = -1;
    uint64_t prev_k_Start_G_addr2 = -1;
    uint64_t prev_k_End_G_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_G_addr1 != -1) {
            if ( calAddrG_addr1( i_Start - prev_i_Start_G_addr1 + prev_i_End_G_addr1, j_Start - prev_j_Start_G_addr1 + prev_j_End_G_addr1, k_Start - prev_k_Start_G_addr1 + prev_k_End_G_addr1) == calAddrG_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_G_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_G_addr2 != -1) {
            if ( calAddrG_addr2( i_Start - prev_i_Start_G_addr2 + prev_i_End_G_addr2, j_Start - prev_j_Start_G_addr2 + prev_j_End_G_addr2, k_Start - prev_k_Start_G_addr2 + prev_k_End_G_addr2) == calAddrG_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_G_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrG_addr0( i, j) == calAddrG_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB8 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr1( i, j, k) == calAddrG_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_G_addr1 = cnt;
                            prev_i_Start_G_addr1 = i_Start;
                            prev_i_End_G_addr1 = i;
                            prev_j_Start_G_addr1 = j_Start;
                            prev_j_End_G_addr1 = j;
                            prev_k_Start_G_addr1 = k_Start;
                            prev_k_End_G_addr1 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr2( i, j, k) == calAddrG_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_G_addr2 = cnt;
                            prev_i_Start_G_addr2 = i_Start;
                            prev_i_End_G_addr2 = i;
                            prev_j_Start_G_addr2 = j_Start;
                            prev_j_End_G_addr2 = j;
                            prev_k_Start_G_addr2 = k_Start;
                            prev_k_End_G_addr2 = k;
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
void ref_G_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_G_addr1 = -1;
    uint64_t prev_i_Start_G_addr1 = -1;
    uint64_t prev_i_End_G_addr1 = -1;
    uint64_t prev_j_Start_G_addr1 = -1;
    uint64_t prev_j_End_G_addr1 = -1;
    uint64_t prev_k_Start_G_addr1 = -1;
    uint64_t prev_k_End_G_addr1 = -1;
    uint64_t prev_cnt_G_addr2 = -1;
    uint64_t prev_i_Start_G_addr2 = -1;
    uint64_t prev_i_End_G_addr2 = -1;
    uint64_t prev_j_Start_G_addr2 = -1;
    uint64_t prev_j_End_G_addr2 = -1;
    uint64_t prev_k_Start_G_addr2 = -1;
    uint64_t prev_k_End_G_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_G_addr1 != -1) {
            if ( calAddrG_addr1( i_Start - prev_i_Start_G_addr1 + prev_i_End_G_addr1, j_Start - prev_j_Start_G_addr1 + prev_j_End_G_addr1, k_Start - prev_k_Start_G_addr1 + prev_k_End_G_addr1) == calAddrG_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_G_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_G_addr2 != -1) {
            if ( calAddrG_addr2( i_Start - prev_i_Start_G_addr2 + prev_i_End_G_addr2, j_Start - prev_j_Start_G_addr2 + prev_j_End_G_addr2, k_Start - prev_k_Start_G_addr2 + prev_k_End_G_addr2) == calAddrG_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_G_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB6 = i_Start;
        for ( int i = iLB6; i < 256; i++) {
            {
            int jLB7 = 0;
            if ( i == i_Start ) {
                jLB7 = j_Start;
            }
            for ( int j = jLB7; j < 256; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrG_addr0( i, j) == calAddrG_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB8 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB8 = k_Start;
                }
                for ( int k = kLB8; k < 256; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr1( i, j, k) == calAddrG_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_G_addr1 = cnt;
                            prev_i_Start_G_addr1 = i_Start;
                            prev_i_End_G_addr1 = i;
                            prev_j_Start_G_addr1 = j_Start;
                            prev_j_End_G_addr1 = j;
                            prev_k_Start_G_addr1 = k_Start;
                            prev_k_End_G_addr1 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr2( i, j, k) == calAddrG_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_G_addr2 = cnt;
                            prev_i_Start_G_addr2 = i_Start;
                            prev_i_End_G_addr2 = i;
                            prev_j_Start_G_addr2 = j_Start;
                            prev_j_End_G_addr2 = j;
                            prev_k_Start_G_addr2 = k_Start;
                            prev_k_End_G_addr2 = k;
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
    ref_A_addr0();
    ref_B_addr0();
    ref_C_addr0();
    ref_D_addr0();
    ref_E_addr0();
    ref_E_addr1();
    ref_E_addr2();
    ref_E_addr3();
    ref_F_addr0();
    ref_F_addr1();
    ref_F_addr2();
    ref_F_addr3();
    ref_G_addr0();
    ref_G_addr1();
    ref_G_addr2();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
mm2 */ 
