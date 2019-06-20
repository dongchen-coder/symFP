
 /* Start to analysis array index
Array index info
mean.addr j
data.addr ((i * 1024) + j)
mean.addr j
mean.addr j
mean.addr j
mean.addr j
mean.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
symmat.addr ((j1 * 1024) + j2)
data.addr ((i * 1024) + j1)
data.addr ((i * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j2 * 1024) + j1)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %m
i32 %n
double* %data
double* %symmat
double* %mean

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--j
--Loop Bound: (0, 1024)
--Loop inc: (j + 1)
--Loop predicate: <
----array access mean.addr j
----i
----Loop Bound: (0, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------array access data.addr ((i * 1024) + j)
------array access mean.addr j
------array access mean.addr j
----array access mean.addr j
----array access mean.addr j
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access mean.addr j
------array access data.addr ((i * 1024) + j)
------array access data.addr ((i * 1024) + j)
--j1
--Loop Bound: (0, 1024)
--Loop inc: (j1 + 1)
--Loop predicate: <
----j2
----Loop Bound: (j1, 1024)
----Loop inc: (j2 + 1)
----Loop predicate: <
------array access symmat.addr ((j1 * 1024) + j2)
------i
------Loop Bound: (0, 1024)
------Loop inc: (i + 1)
------Loop predicate: <
--------array access data.addr ((i * 1024) + j1)
--------array access data.addr ((i * 1024) + j2)
--------array access symmat.addr ((j1 * 1024) + j2)
--------array access symmat.addr ((j1 * 1024) + j2)
------array access symmat.addr ((j1 * 1024) + j2)
------array access symmat.addr ((j2 * 1024) + j1)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 0 0 
Dump stride: 1 1 
init counter: 0 0 0 
Dump stride: 1 1 1 
Dump tree:
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 52
--------Sample number: 537
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
int calAddrmean_addr0( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrdata_addr0( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrmean_addr1( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrmean_addr2( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrmean_addr3( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrmean_addr4( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrmean_addr5( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrdata_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrdata_addr2( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr0( int j1, int j2) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrdata_addr3( int j1, int j2, int i) {
    int result = (((i * 1024) + j1)) * 8 / 64;
    return result;
}
int calAddrdata_addr4( int j1, int j2, int i) {
    int result = (((i * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr1( int j1, int j2, int i) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr2( int j1, int j2, int i) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr3( int j1, int j2) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr4( int j1, int j2) {
    int result = (((j2 * 1024) + j1)) * 8 / 64;
    return result;
}
void ref_data_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB0 = j_Start;
        for ( int j = jLB0; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB1 = 0;
            if ( j == j_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr0( j, i) == calAddrdata_addr0(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
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
                    if ( calAddrdata_addr1( i, j) == calAddrdata_addr0(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr2( i, j) == calAddrdata_addr0(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int j1LB4 = 0;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr3( j1, j2, i) == calAddrdata_addr0(j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr4( j1, j2, i) == calAddrdata_addr0(j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
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
void ref_data_addr1() {
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
                    if ( calAddrdata_addr1( i, j) == calAddrdata_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr2( i, j) == calAddrdata_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int j1LB4 = 0;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr3( j1, j2, i) == calAddrdata_addr1(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr4( j1, j2, i) == calAddrdata_addr1(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
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
void ref_data_addr2() {
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
                    if ( calAddrdata_addr1( i, j) == calAddrdata_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr2( i, j) == calAddrdata_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        {
        int j1LB4 = 0;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr3( j1, j2, i) == calAddrdata_addr2(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr4( j1, j2, i) == calAddrdata_addr2(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
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
void ref_data_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 537;) {
SAMPLE:
        int j1_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - j1_Start) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - j1_Start) + j1_Start;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB4 = j1_Start;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            if ( j1 == j1_Start ) {
                j2LB5 = j2_Start;
            }
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                if ( j1 == j1_Start && j2 == j2_Start ) {
                    iLB6 = i_Start;
                }
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr3( j1, j2, i) == calAddrdata_addr3(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr4( j1, j2, i) == calAddrdata_addr3(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
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
void ref_data_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 537;) {
SAMPLE:
        int j1_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - j1_Start) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - j1_Start) + j1_Start;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB4 = j1_Start;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            if ( j1 == j1_Start ) {
                j2LB5 = j2_Start;
            }
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                if ( j1 == j1_Start && j2 == j2_Start ) {
                    iLB6 = i_Start;
                }
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr3( j1, j2, i) == calAddrdata_addr4(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr4( j1, j2, i) == calAddrdata_addr4(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
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
void ref_mean_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB0 = j_Start;
        for ( int j = jLB0; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr0( j) == calAddrmean_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr1( j, i) == calAddrmean_addr0(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr0(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr3( j) == calAddrmean_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( i, j) == calAddrmean_addr0(j_Start)) {
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
        int j1LB4 = 0;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
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
void ref_mean_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB0 = j_Start;
        for ( int j = jLB0; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr0( j) == calAddrmean_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB1 = 0;
            if ( j == j_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr1( j, i) == calAddrmean_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr3( j) == calAddrmean_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( i, j) == calAddrmean_addr1(j_Start, i_Start)) {
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
        int j1LB4 = 0;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
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
void ref_mean_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB0 = j_Start;
        for ( int j = jLB0; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr0( j) == calAddrmean_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB1 = 0;
            if ( j == j_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr1( j, i) == calAddrmean_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr3( j) == calAddrmean_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( i, j) == calAddrmean_addr2(j_Start, i_Start)) {
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
        int j1LB4 = 0;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
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
void ref_mean_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB0 = j_Start;
        for ( int j = jLB0; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr0( j) == calAddrmean_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr1( j, i) == calAddrmean_addr3(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr3(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr3( j) == calAddrmean_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( i, j) == calAddrmean_addr3(j_Start)) {
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
        int j1LB4 = 0;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
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
void ref_mean_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB0 = j_Start;
        for ( int j = jLB0; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr0( j) == calAddrmean_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr1( j, i) == calAddrmean_addr4(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr4(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr3( j) == calAddrmean_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
        }
        }
        {
        int iLB2 = 0;
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 0;
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( i, j) == calAddrmean_addr4(j_Start)) {
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
        int j1LB4 = 0;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
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
void ref_mean_addr5() {
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
                    if ( calAddrmean_addr5( i, j) == calAddrmean_addr5(i_Start, j_Start)) {
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
        {
        int j1LB4 = 0;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
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
void ref_symmat_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int j1_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - j1_Start) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - j1_Start) + j1_Start;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB4 = j1_Start;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            if ( j1 == j1_Start ) {
                j2LB5 = j2_Start;
            }
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr0( j1, j2) == calAddrsymmat_addr0(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr1( j1, j2, i) == calAddrsymmat_addr0(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr0(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr3( j1, j2) == calAddrsymmat_addr0(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr0(j1_Start, j2_Start)) {
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
void ref_symmat_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 537;) {
SAMPLE:
        int j1_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - j1_Start) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - j1_Start) + j1_Start;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB4 = j1_Start;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            if ( j1 == j1_Start ) {
                j2LB5 = j2_Start;
            }
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr0( j1, j2) == calAddrsymmat_addr1(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB6 = 0;
                if ( j1 == j1_Start && j2 == j2_Start ) {
                    iLB6 = i_Start;
                }
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr1( j1, j2, i) == calAddrsymmat_addr1(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr1(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr3( j1, j2) == calAddrsymmat_addr1(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr1(j1_Start, j2_Start, i_Start)) {
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
void ref_symmat_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 537;) {
SAMPLE:
        int j1_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - j1_Start) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - j1_Start) + j1_Start;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB4 = j1_Start;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            if ( j1 == j1_Start ) {
                j2LB5 = j2_Start;
            }
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr0( j1, j2) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB6 = 0;
                if ( j1 == j1_Start && j2 == j2_Start ) {
                    iLB6 = i_Start;
                }
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr1( j1, j2, i) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr3( j1, j2) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
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
void ref_symmat_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int j1_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - j1_Start) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - j1_Start) + j1_Start;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB4 = j1_Start;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            if ( j1 == j1_Start ) {
                j2LB5 = j2_Start;
            }
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr0( j1, j2) == calAddrsymmat_addr3(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr1( j1, j2, i) == calAddrsymmat_addr3(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr3(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr3( j1, j2) == calAddrsymmat_addr3(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr3(j1_Start, j2_Start)) {
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
void ref_symmat_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int j1_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - j1_Start) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - j1_Start) + j1_Start;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB4 = j1_Start;
        for ( int j1 = j1LB4; j1 < 1024; j1++) {
            {
            int j2LB5 = j1;
            if ( j1 == j1_Start ) {
                j2LB5 = j2_Start;
            }
            for ( int j2 = j2LB5; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr0( j1, j2) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB6 = 0;
                for ( int i = iLB6; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr1( j1, j2, i) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr3( j1, j2) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
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
    ref_data_addr0();
    ref_data_addr1();
    ref_data_addr2();
    ref_data_addr3();
    ref_data_addr4();
    ref_mean_addr0();
    ref_mean_addr1();
    ref_mean_addr2();
    ref_mean_addr3();
    ref_mean_addr4();
    ref_mean_addr5();
    ref_symmat_addr0();
    ref_symmat_addr1();
    ref_symmat_addr2();
    ref_symmat_addr3();
    ref_symmat_addr4();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: covariance */ 
