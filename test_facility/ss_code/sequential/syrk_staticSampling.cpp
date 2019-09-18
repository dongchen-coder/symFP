
 /* Start to analysis array index
Array index info: Total number of references: 6
C.addr ((i * 256) + l)
C.addr ((i * 256) + l)
C.addr ((i * 256) + j)
C.addr ((i * 256) + j)
A.addr ((i * 256) + l)
A.addr ((j * 256) + l)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %ni
i32 %nj
double %alpha
double %beta
double* %A
double* %C

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 256)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, i)
----Loop inc: (j + 1)
----Loop predicate: <=
------array access C.addr ((i * 256) + j)
------array access C.addr ((i * 256) + j)
------k
------Loop Bound: (0, 256)
------Loop inc: (k + 1)
------Loop predicate: <
--------l
--------Loop Bound: (0, i)
--------Loop inc: (l + 1)
--------Loop predicate: <=
----------array access A.addr ((i * 256) + l)
----------array access A.addr ((j * 256) + l)
----------array access C.addr ((i * 256) + l)
----------array access C.addr ((i * 256) + l)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 0 0 
Dump stride: 1 1 
init counter: 0 0 0 
Dump stride: 1 1 1 
init counter: 0 0 0 0 
Dump stride: 1 1 1 1 
Dump tree:
----Sample number: 2
------Sample number: 3
--------Sample number: 8
----------Sample number: 14
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
/* C_addr ((i * 256) + j) 0 */
int calAddrC_addr0( int i, int j) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* C_addr ((i * 256) + j) 1 */
int calAddrC_addr1( int i, int j) {
    int result = (((i * 256) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 256) + l) 2 */
int calAddrA_addr2( int i, int j, int k, int l) {
    int result = (((i * 256) + l)) * 8 / 64;
    return result;
}
/* A_addr ((j * 256) + l) 3 */
int calAddrA_addr3( int i, int j, int k, int l) {
    int result = (((j * 256) + l)) * 8 / 64;
    return result;
}
/* C_addr ((i * 256) + l) 4 */
int calAddrC_addr4( int i, int j, int k, int l) {
    int result = (((i * 256) + l)) * 8 / 64;
    return result;
}
/* C_addr ((i * 256) + l) 5 */
int calAddrC_addr5( int i, int j, int k, int l) {
    int result = (((i * 256) + l)) * 8 / 64;
    return result;
}
void ref_C_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 14;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0 + 1) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int l_Start = rand() % (i_Start - 0 + 1) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" + std::to_string(l_Start) + "_" ;
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
            for ( int j = jLB1; j <= i; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrC_addr0( i, j) == calAddrC_addr4(i_Start, j_Start, k_Start, l_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrC_addr1( i, j) == calAddrC_addr4(i_Start, j_Start, k_Start, l_Start)) {
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
                    {
                    int lLB3 = 0;
                    if ( i == i_Start && j == j_Start && k == k_Start ) {
                        lLB3 = l_Start;
                    }
                    for ( int l = lLB3; l <= i; l++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC_addr4( i, j, k, l) == calAddrC_addr4(i_Start, j_Start, k_Start, l_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC_addr5( i, j, k, l) == calAddrC_addr4(i_Start, j_Start, k_Start, l_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_C_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 14;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0 + 1) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int l_Start = rand() % (i_Start - 0 + 1) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" + std::to_string(l_Start) + "_" ;
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
            for ( int j = jLB1; j <= i; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrC_addr0( i, j) == calAddrC_addr5(i_Start, j_Start, k_Start, l_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrC_addr1( i, j) == calAddrC_addr5(i_Start, j_Start, k_Start, l_Start)) {
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
                    {
                    int lLB3 = 0;
                    if ( i == i_Start && j == j_Start && k == k_Start ) {
                        lLB3 = l_Start;
                    }
                    for ( int l = lLB3; l <= i; l++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC_addr4( i, j, k, l) == calAddrC_addr5(i_Start, j_Start, k_Start, l_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC_addr5( i, j, k, l) == calAddrC_addr5(i_Start, j_Start, k_Start, l_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_C_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 3;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0 + 1) + 0;
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
            for ( int j = jLB1; j <= i; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrC_addr0( i, j) == calAddrC_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrC_addr1( i, j) == calAddrC_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < 256; k++) {
                    {
                    int lLB3 = 0;
                    for ( int l = lLB3; l <= i; l++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC_addr4( i, j, k, l) == calAddrC_addr0(i_Start, j_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC_addr5( i, j, k, l) == calAddrC_addr0(i_Start, j_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_C_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 3;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0 + 1) + 0;
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
            for ( int j = jLB1; j <= i; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrC_addr0( i, j) == calAddrC_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrC_addr1( i, j) == calAddrC_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < 256; k++) {
                    {
                    int lLB3 = 0;
                    for ( int l = lLB3; l <= i; l++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC_addr4( i, j, k, l) == calAddrC_addr1(i_Start, j_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC_addr5( i, j, k, l) == calAddrC_addr1(i_Start, j_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 14;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0 + 1) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int l_Start = rand() % (i_Start - 0 + 1) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" + std::to_string(l_Start) + "_" ;
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
            for ( int j = jLB1; j <= i; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 256; k++) {
                    {
                    int lLB3 = 0;
                    if ( i == i_Start && j == j_Start && k == k_Start ) {
                        lLB3 = l_Start;
                    }
                    for ( int l = lLB3; l <= i; l++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( i, j, k, l) == calAddrA_addr2(i_Start, j_Start, k_Start, l_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( i, j, k, l) == calAddrA_addr2(i_Start, j_Start, k_Start, l_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 14;) {
SAMPLE:
        int i_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0 + 1) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (256 - 0) + 0;
        if ( (i_Start - 0 + 1) == 0) goto SAMPLE;
        int l_Start = rand() % (i_Start - 0 + 1) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" + std::to_string(l_Start) + "_" ;
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
            for ( int j = jLB1; j <= i; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 256; k++) {
                    {
                    int lLB3 = 0;
                    if ( i == i_Start && j == j_Start && k == k_Start ) {
                        lLB3 = l_Start;
                    }
                    for ( int l = lLB3; l <= i; l++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( i, j, k, l) == calAddrA_addr3(i_Start, j_Start, k_Start, l_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( i, j, k, l) == calAddrA_addr3(i_Start, j_Start, k_Start, l_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
int main() {
    ref_C_addr4();
    ref_C_addr5();
    ref_C_addr0();
    ref_C_addr1();
    ref_A_addr2();
    ref_A_addr3();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: syrk */ 
