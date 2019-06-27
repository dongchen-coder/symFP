
 /* Start to analysis array index
Array index info
A.addr ((i * 1024) + j)
u1.addr i
v1.addr j
u2.addr i
v2.addr j
A.addr ((i * 1024) + j)
x.addr i
A.addr ((j * 1024) + i)
y.addr j
x.addr i
x.addr i
z.addr i
x.addr i
w.addr i
A.addr ((i * 1024) + j)
x.addr j
w.addr i

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
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access A.addr ((i * 1024) + j)
------array access u1.addr i
------array access v1.addr j
------array access u2.addr i
------array access v2.addr j
------array access A.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access x.addr i
------array access A.addr ((j * 1024) + i)
------array access y.addr j
------array access x.addr i
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----array access x.addr i
----array access z.addr i
----array access x.addr i
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access w.addr i
------array access A.addr ((i * 1024) + j)
------array access x.addr j
------array access w.addr i

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 1
------Sample number: 1
----Sample number: 1
------Sample number: 1
----Sample number: 1
----Sample number: 1
------Sample number: 1
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
int calAddrA_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddru1_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrv1_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddru2_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrv2_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrA_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrx_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrA_addr2( int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 64;
    return result;
}
int calAddry_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrx_addr1( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrx_addr2( int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrz_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrx_addr3( int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrw_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrA_addr3( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrx_addr4( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrw_addr1( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0) == calAddrA_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1) == calAddrA_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr0 = cnt;
                        prev_i_Start_A_addr0 = i_Start;
                        prev_i_End_A_addr0 = i;
                        prev_j_Start_A_addr0 = j_Start;
                        prev_j_End_A_addr0 = j;
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
                    if ( calAddrA_addr1( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr1 = cnt;
                        prev_i_Start_A_addr1 = i_Start;
                        prev_i_End_A_addr1 = i;
                        prev_j_Start_A_addr1 = j_Start;
                        prev_j_End_A_addr1 = j;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB2 = 0;
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( i, j) == calAddrA_addr0(i_Start, j_Start)) {
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
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( i, j) == calAddrA_addr0(i_Start, j_Start)) {
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
EndSample:
        s++;
        }
}
void ref_A_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0) == calAddrA_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1) == calAddrA_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( i, j) == calAddrA_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr0 = cnt;
                        prev_i_Start_A_addr0 = i_Start;
                        prev_i_End_A_addr0 = i;
                        prev_j_Start_A_addr0 = j_Start;
                        prev_j_End_A_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( i, j) == calAddrA_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr1 = cnt;
                        prev_i_Start_A_addr1 = i_Start;
                        prev_i_End_A_addr1 = i;
                        prev_j_Start_A_addr1 = j_Start;
                        prev_j_End_A_addr1 = j;
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
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( i, j) == calAddrA_addr1(i_Start, j_Start)) {
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
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( i, j) == calAddrA_addr1(i_Start, j_Start)) {
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
EndSample:
        s++;
        }
}
void ref_A_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2) == calAddrA_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( i, j) == calAddrA_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr2 = cnt;
                        prev_i_Start_A_addr2 = i_Start;
                        prev_i_End_A_addr2 = i;
                        prev_j_Start_A_addr2 = j_Start;
                        prev_j_End_A_addr2 = j;
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
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( i, j) == calAddrA_addr2(i_Start, j_Start)) {
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
EndSample:
        s++;
        }
}
void ref_A_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3) == calAddrA_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( i, j) == calAddrA_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr3 = cnt;
                        prev_i_Start_A_addr3 = i_Start;
                        prev_i_End_A_addr3 = i;
                        prev_j_Start_A_addr3 = j_Start;
                        prev_j_End_A_addr3 = j;
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
void ref_u1_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u1_addr0 = -1;
    uint64_t prev_i_Start_u1_addr0 = -1;
    uint64_t prev_i_End_u1_addr0 = -1;
    uint64_t prev_j_Start_u1_addr0 = -1;
    uint64_t prev_j_End_u1_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_u1_addr0 != -1) {
            if ( calAddru1_addr0( i_Start - prev_i_Start_u1_addr0 + prev_i_End_u1_addr0, j_Start - prev_j_Start_u1_addr0 + prev_j_End_u1_addr0) == calAddru1_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u1_addr0);
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
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru1_addr0( i, j) == calAddru1_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u1_addr0 = cnt;
                        prev_i_Start_u1_addr0 = i_Start;
                        prev_i_End_u1_addr0 = i;
                        prev_j_Start_u1_addr0 = j_Start;
                        prev_j_End_u1_addr0 = j;
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
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 1024; j++) {
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
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
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
void ref_u2_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u2_addr0 = -1;
    uint64_t prev_i_Start_u2_addr0 = -1;
    uint64_t prev_i_End_u2_addr0 = -1;
    uint64_t prev_j_Start_u2_addr0 = -1;
    uint64_t prev_j_End_u2_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_u2_addr0 != -1) {
            if ( calAddru2_addr0( i_Start - prev_i_Start_u2_addr0 + prev_i_End_u2_addr0, j_Start - prev_j_Start_u2_addr0 + prev_j_End_u2_addr0) == calAddru2_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u2_addr0);
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
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru2_addr0( i, j) == calAddru2_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u2_addr0 = cnt;
                        prev_i_Start_u2_addr0 = i_Start;
                        prev_i_End_u2_addr0 = i;
                        prev_j_Start_u2_addr0 = j_Start;
                        prev_j_End_u2_addr0 = j;
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
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 1024; j++) {
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
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
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
void ref_v1_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v1_addr0 = -1;
    uint64_t prev_i_Start_v1_addr0 = -1;
    uint64_t prev_i_End_v1_addr0 = -1;
    uint64_t prev_j_Start_v1_addr0 = -1;
    uint64_t prev_j_End_v1_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_v1_addr0 != -1) {
            if ( calAddrv1_addr0( i_Start - prev_i_Start_v1_addr0 + prev_i_End_v1_addr0, j_Start - prev_j_Start_v1_addr0 + prev_j_End_v1_addr0) == calAddrv1_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v1_addr0);
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
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv1_addr0( i, j) == calAddrv1_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v1_addr0 = cnt;
                        prev_i_Start_v1_addr0 = i_Start;
                        prev_i_End_v1_addr0 = i;
                        prev_j_Start_v1_addr0 = j_Start;
                        prev_j_End_v1_addr0 = j;
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
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 1024; j++) {
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
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
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
void ref_v2_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v2_addr0 = -1;
    uint64_t prev_i_Start_v2_addr0 = -1;
    uint64_t prev_i_End_v2_addr0 = -1;
    uint64_t prev_j_Start_v2_addr0 = -1;
    uint64_t prev_j_End_v2_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_v2_addr0 != -1) {
            if ( calAddrv2_addr0( i_Start - prev_i_Start_v2_addr0 + prev_i_End_v2_addr0, j_Start - prev_j_Start_v2_addr0 + prev_j_End_v2_addr0) == calAddrv2_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v2_addr0);
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
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv2_addr0( i, j) == calAddrv2_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v2_addr0 = cnt;
                        prev_i_Start_v2_addr0 = i_Start;
                        prev_i_End_v2_addr0 = i;
                        prev_j_Start_v2_addr0 = j_Start;
                        prev_j_End_v2_addr0 = j;
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
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 1024; j++) {
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
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
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
void ref_w_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_w_addr0 = -1;
    uint64_t prev_i_Start_w_addr0 = -1;
    uint64_t prev_i_End_w_addr0 = -1;
    uint64_t prev_j_Start_w_addr0 = -1;
    uint64_t prev_j_End_w_addr0 = -1;
    uint64_t prev_cnt_w_addr1 = -1;
    uint64_t prev_i_Start_w_addr1 = -1;
    uint64_t prev_i_End_w_addr1 = -1;
    uint64_t prev_j_Start_w_addr1 = -1;
    uint64_t prev_j_End_w_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_w_addr0 != -1) {
            if ( calAddrw_addr0( i_Start - prev_i_Start_w_addr0 + prev_i_End_w_addr0, j_Start - prev_j_Start_w_addr0 + prev_j_End_w_addr0) == calAddrw_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_w_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_w_addr1 != -1) {
            if ( calAddrw_addr1( i_Start - prev_i_Start_w_addr1 + prev_i_End_w_addr1, j_Start - prev_j_Start_w_addr1 + prev_j_End_w_addr1) == calAddrw_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_w_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr0( i, j) == calAddrw_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_w_addr0 = cnt;
                        prev_i_Start_w_addr0 = i_Start;
                        prev_i_End_w_addr0 = i;
                        prev_j_Start_w_addr0 = j_Start;
                        prev_j_End_w_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr1( i, j) == calAddrw_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_w_addr1 = cnt;
                        prev_i_Start_w_addr1 = i_Start;
                        prev_i_End_w_addr1 = i;
                        prev_j_Start_w_addr1 = j_Start;
                        prev_j_End_w_addr1 = j;
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
void ref_w_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_w_addr0 = -1;
    uint64_t prev_i_Start_w_addr0 = -1;
    uint64_t prev_i_End_w_addr0 = -1;
    uint64_t prev_j_Start_w_addr0 = -1;
    uint64_t prev_j_End_w_addr0 = -1;
    uint64_t prev_cnt_w_addr1 = -1;
    uint64_t prev_i_Start_w_addr1 = -1;
    uint64_t prev_i_End_w_addr1 = -1;
    uint64_t prev_j_Start_w_addr1 = -1;
    uint64_t prev_j_End_w_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_w_addr0 != -1) {
            if ( calAddrw_addr0( i_Start - prev_i_Start_w_addr0 + prev_i_End_w_addr0, j_Start - prev_j_Start_w_addr0 + prev_j_End_w_addr0) == calAddrw_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_w_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_w_addr1 != -1) {
            if ( calAddrw_addr1( i_Start - prev_i_Start_w_addr1 + prev_i_End_w_addr1, j_Start - prev_j_Start_w_addr1 + prev_j_End_w_addr1) == calAddrw_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_w_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr0( i, j) == calAddrw_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_w_addr0 = cnt;
                        prev_i_Start_w_addr0 = i_Start;
                        prev_i_End_w_addr0 = i;
                        prev_j_Start_w_addr0 = j_Start;
                        prev_j_End_w_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr1( i, j) == calAddrw_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_w_addr1 = cnt;
                        prev_i_Start_w_addr1 = i_Start;
                        prev_i_End_w_addr1 = i;
                        prev_j_Start_w_addr1 = j_Start;
                        prev_j_End_w_addr1 = j;
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
void ref_x_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr0 = -1;
    uint64_t prev_i_Start_x_addr0 = -1;
    uint64_t prev_i_End_x_addr0 = -1;
    uint64_t prev_j_Start_x_addr0 = -1;
    uint64_t prev_j_End_x_addr0 = -1;
    uint64_t prev_cnt_x_addr1 = -1;
    uint64_t prev_i_Start_x_addr1 = -1;
    uint64_t prev_i_End_x_addr1 = -1;
    uint64_t prev_j_Start_x_addr1 = -1;
    uint64_t prev_j_End_x_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_x_addr0 != -1) {
            if ( calAddrx_addr0( i_Start - prev_i_Start_x_addr0 + prev_i_End_x_addr0, j_Start - prev_j_Start_x_addr0 + prev_j_End_x_addr0) == calAddrx_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr1 != -1) {
            if ( calAddrx_addr1( i_Start - prev_i_Start_x_addr1 + prev_i_End_x_addr1, j_Start - prev_j_Start_x_addr1 + prev_j_End_x_addr1) == calAddrx_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr0( i, j) == calAddrx_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr0 = cnt;
                        prev_i_Start_x_addr0 = i_Start;
                        prev_i_End_x_addr0 = i;
                        prev_j_Start_x_addr0 = j_Start;
                        prev_j_End_x_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr1 = cnt;
                        prev_i_Start_x_addr1 = i_Start;
                        prev_i_End_x_addr1 = i;
                        prev_j_Start_x_addr1 = j_Start;
                        prev_j_End_x_addr1 = j;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( i) == calAddrx_addr0(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( i) == calAddrx_addr0(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( i, j) == calAddrx_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
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
void ref_x_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr0 = -1;
    uint64_t prev_i_Start_x_addr0 = -1;
    uint64_t prev_i_End_x_addr0 = -1;
    uint64_t prev_j_Start_x_addr0 = -1;
    uint64_t prev_j_End_x_addr0 = -1;
    uint64_t prev_cnt_x_addr1 = -1;
    uint64_t prev_i_Start_x_addr1 = -1;
    uint64_t prev_i_End_x_addr1 = -1;
    uint64_t prev_j_Start_x_addr1 = -1;
    uint64_t prev_j_End_x_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_x_addr0 != -1) {
            if ( calAddrx_addr0( i_Start - prev_i_Start_x_addr0 + prev_i_End_x_addr0, j_Start - prev_j_Start_x_addr0 + prev_j_End_x_addr0) == calAddrx_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr1 != -1) {
            if ( calAddrx_addr1( i_Start - prev_i_Start_x_addr1 + prev_i_End_x_addr1, j_Start - prev_j_Start_x_addr1 + prev_j_End_x_addr1) == calAddrx_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr0( i, j) == calAddrx_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr0 = cnt;
                        prev_i_Start_x_addr0 = i_Start;
                        prev_i_End_x_addr0 = i;
                        prev_j_Start_x_addr0 = j_Start;
                        prev_j_End_x_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr1 = cnt;
                        prev_i_Start_x_addr1 = i_Start;
                        prev_i_End_x_addr1 = i;
                        prev_j_Start_x_addr1 = j_Start;
                        prev_j_End_x_addr1 = j;
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
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( i) == calAddrx_addr1(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( i) == calAddrx_addr1(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( i, j) == calAddrx_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
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
void ref_x_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr2 = -1;
    uint64_t prev_i_Start_x_addr2 = -1;
    uint64_t prev_i_End_x_addr2 = -1;
    uint64_t prev_cnt_x_addr3 = -1;
    uint64_t prev_i_Start_x_addr3 = -1;
    uint64_t prev_i_End_x_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_x_addr2 != -1) {
            if ( calAddrx_addr2( i_Start - prev_i_Start_x_addr2 + prev_i_End_x_addr2) == calAddrx_addr2(i_Start)) {
                rtHistoCal(prev_cnt_x_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr3 != -1) {
            if ( calAddrx_addr3( i_Start - prev_i_Start_x_addr3 + prev_i_End_x_addr3) == calAddrx_addr2(i_Start)) {
                rtHistoCal(prev_cnt_x_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( i) == calAddrx_addr2(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr2 = cnt;
                    prev_i_Start_x_addr2 = i_Start;
                    prev_i_End_x_addr2 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( i) == calAddrx_addr2(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr3 = cnt;
                    prev_i_Start_x_addr3 = i_Start;
                    prev_i_End_x_addr3 = i;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( i, j) == calAddrx_addr2(i_Start)) {
                        rtHistoCal(cnt);
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
void ref_x_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr2 = -1;
    uint64_t prev_i_Start_x_addr2 = -1;
    uint64_t prev_i_End_x_addr2 = -1;
    uint64_t prev_cnt_x_addr3 = -1;
    uint64_t prev_i_Start_x_addr3 = -1;
    uint64_t prev_i_End_x_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_x_addr2 != -1) {
            if ( calAddrx_addr2( i_Start - prev_i_Start_x_addr2 + prev_i_End_x_addr2) == calAddrx_addr3(i_Start)) {
                rtHistoCal(prev_cnt_x_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr3 != -1) {
            if ( calAddrx_addr3( i_Start - prev_i_Start_x_addr3 + prev_i_End_x_addr3) == calAddrx_addr3(i_Start)) {
                rtHistoCal(prev_cnt_x_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( i) == calAddrx_addr3(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr2 = cnt;
                    prev_i_Start_x_addr2 = i_Start;
                    prev_i_End_x_addr2 = i;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( i) == calAddrx_addr3(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr3 = cnt;
                    prev_i_Start_x_addr3 = i_Start;
                    prev_i_End_x_addr3 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( i, j) == calAddrx_addr3(i_Start)) {
                        rtHistoCal(cnt);
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
void ref_x_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr4 = -1;
    uint64_t prev_i_Start_x_addr4 = -1;
    uint64_t prev_i_End_x_addr4 = -1;
    uint64_t prev_j_Start_x_addr4 = -1;
    uint64_t prev_j_End_x_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_x_addr4 != -1) {
            if ( calAddrx_addr4( i_Start - prev_i_Start_x_addr4 + prev_i_End_x_addr4, j_Start - prev_j_Start_x_addr4 + prev_j_End_x_addr4) == calAddrx_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( i, j) == calAddrx_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr4 = cnt;
                        prev_i_Start_x_addr4 = i_Start;
                        prev_i_End_x_addr4 = i;
                        prev_j_Start_x_addr4 = j_Start;
                        prev_j_End_x_addr4 = j;
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
void ref_y_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_y_addr0 = -1;
    uint64_t prev_i_Start_y_addr0 = -1;
    uint64_t prev_i_End_y_addr0 = -1;
    uint64_t prev_j_Start_y_addr0 = -1;
    uint64_t prev_j_End_y_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
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
        if ( prev_cnt_y_addr0 != -1) {
            if ( calAddry_addr0( i_Start - prev_i_Start_y_addr0 + prev_i_End_y_addr0, j_Start - prev_j_Start_y_addr0 + prev_j_End_y_addr0) == calAddry_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_y_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr0( i, j) == calAddry_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_y_addr0 = cnt;
                        prev_i_Start_y_addr0 = i_Start;
                        prev_i_End_y_addr0 = i;
                        prev_j_Start_y_addr0 = j_Start;
                        prev_j_End_y_addr0 = j;
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
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
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
void ref_z_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_z_addr0 = -1;
    uint64_t prev_i_Start_z_addr0 = -1;
    uint64_t prev_i_End_z_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_z_addr0 != -1) {
            if ( calAddrz_addr0( i_Start - prev_i_Start_z_addr0 + prev_i_End_z_addr0) == calAddrz_addr0(i_Start)) {
                rtHistoCal(prev_cnt_z_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrz_addr0( i) == calAddrz_addr0(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_z_addr0 = cnt;
                    prev_i_Start_z_addr0 = i_Start;
                    prev_i_End_z_addr0 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB5 = 0;
        for ( int i = iLB5; i < 1024; i++) {
            {
            int jLB6 = 0;
            for ( int j = jLB6; j < 1024; j++) {
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
int main() {
    ref_A_addr0();
    ref_A_addr1();
    ref_A_addr2();
    ref_A_addr3();
    ref_u1_addr0();
    ref_u2_addr0();
    ref_v1_addr0();
    ref_v2_addr0();
    ref_w_addr0();
    ref_w_addr1();
    ref_x_addr0();
    ref_x_addr1();
    ref_x_addr2();
    ref_x_addr3();
    ref_x_addr4();
    ref_y_addr0();
    ref_z_addr0();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: gemver */ 
