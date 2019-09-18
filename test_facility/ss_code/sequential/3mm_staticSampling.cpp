
 /* Start to analysis array index
Array index info: Total number of references: 15
E.addr ((i * 256) + j)
E.addr ((i * 256) + j)
E.addr ((i * 256) + j)
C.addr ((i * 256) + k)
D.addr ((k * 256) + j)
A.addr ((i * 256) + k)
B.addr ((k * 256) + j)
G.addr ((i * 256) + j)
E.addr ((i * 256) + k)
F.addr ((k * 256) + j)
F.addr ((i * 256) + j)
F.addr ((i * 256) + j)
F.addr ((i * 256) + j)
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
/* E_addr ((i * 256) + j) 0 */
int calAddrE_addr0( int i, int j) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 256) + k) 1 */
int calAddrA_addr1( int i, int j, int k) {
    int result = (((i * 256) + k)) * 8 / 64;
    return result;
}
/* B_addr ((k * 256) + j) 2 */
int calAddrB_addr2( int i, int j, int k) {
    int result = (((k * 256) + j)) * 8 / 64;
    return result;
}
/* E_addr ((i * 256) + j) 3 */
int calAddrE_addr3( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* E_addr ((i * 256) + j) 4 */
int calAddrE_addr4( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* F_addr ((i * 256) + j) 5 */
int calAddrF_addr5( int i, int j) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* C_addr ((i * 256) + k) 6 */
int calAddrC_addr6( int i, int j, int k) {
    int result = (((i * 256) + k)) * 8 / 64;
    return result;
}
/* D_addr ((k * 256) + j) 7 */
int calAddrD_addr7( int i, int j, int k) {
    int result = (((k * 256) + j)) * 8 / 64;
    return result;
}
/* F_addr ((i * 256) + j) 8 */
int calAddrF_addr8( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* F_addr ((i * 256) + j) 9 */
int calAddrF_addr9( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* G_addr ((i * 256) + j) 10 */
int calAddrG_addr10( int i, int j) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* E_addr ((i * 256) + k) 11 */
int calAddrE_addr11( int i, int j, int k) {
    int result = (((i * 256) + k)) * 8 / 64;
    return result;
}
/* F_addr ((k * 256) + j) 12 */
int calAddrF_addr12( int i, int j, int k) {
    int result = (((k * 256) + j)) * 8 / 64;
    return result;
}
/* G_addr ((i * 256) + j) 13 */
int calAddrG_addr13( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* G_addr ((i * 256) + j) 14 */
int calAddrG_addr14( int i, int j, int k) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
void ref_E_addr3() {
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
                    if ( calAddrE_addr0( i, j) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr4( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
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
void ref_E_addr4() {
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
                    if ( calAddrE_addr0( i, j) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr4( i, j, k) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr4(i_Start, j_Start, k_Start)) {
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
void ref_E_addr0() {
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
                        if ( calAddrE_addr3( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrE_addr4( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
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
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr0(i_Start, j_Start)) {
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
void ref_C_addr6() {
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
                        if ( calAddrC_addr6( i, j, k) == calAddrC_addr6(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
void ref_D_addr7() {
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
                        if ( calAddrD_addr7( i, j, k) == calAddrD_addr7(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
void ref_A_addr1() {
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
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
void ref_B_addr2() {
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
                        if ( calAddrB_addr2( i, j, k) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
void ref_G_addr10() {
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
                    if ( calAddrG_addr10( i, j) == calAddrG_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
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
                        if ( calAddrG_addr13( i, j, k) == calAddrG_addr10(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr14( i, j, k) == calAddrG_addr10(i_Start, j_Start)) {
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
void ref_E_addr11() {
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
                        if ( calAddrE_addr11( i, j, k) == calAddrE_addr11(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr12(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
                    if ( calAddrF_addr5( i, j) == calAddrF_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
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
                        if ( calAddrF_addr8( i, j, k) == calAddrF_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr9( i, j, k) == calAddrF_addr5(i_Start, j_Start)) {
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
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr5(i_Start, j_Start)) {
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
void ref_F_addr8() {
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
                    if ( calAddrF_addr5( i, j) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrF_addr8( i, j, k) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr9( i, j, k) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr8(i_Start, j_Start, k_Start)) {
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
void ref_F_addr9() {
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
                    if ( calAddrF_addr5( i, j) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrF_addr8( i, j, k) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrF_addr9( i, j, k) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
                        if ( calAddrF_addr12( i, j, k) == calAddrF_addr9(i_Start, j_Start, k_Start)) {
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
void ref_G_addr13() {
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
                    if ( calAddrG_addr10( i, j) == calAddrG_addr13(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrG_addr13( i, j, k) == calAddrG_addr13(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr14( i, j, k) == calAddrG_addr13(i_Start, j_Start, k_Start)) {
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
void ref_G_addr14() {
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
                    if ( calAddrG_addr10( i, j) == calAddrG_addr14(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrG_addr13( i, j, k) == calAddrG_addr14(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrG_addr14( i, j, k) == calAddrG_addr14(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
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
    return 0;
}
 /* Analyze function: mm2 */ 
