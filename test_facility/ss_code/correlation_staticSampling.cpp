
 /* Start to analysis array index
Array index info
mean.addr j
data.addr ((i * 1024) + j)
mean.addr j
mean.addr j
mean.addr j
mean.addr j
stddev.addr j
data.addr ((i * 1024) + j)
mean.addr j
data.addr ((i * 1024) + j)
mean.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
mean.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
stddev.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
symmat.addr ((j1 * 1024) + j1)
symmat.addr ((j1 * 1024) + j2)
data.addr ((i * 1024) + j1)
data.addr ((i * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j2 * 1024) + j1)
symmat.addr 1048575

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %m
i32 %n
double* %data
double* %mean
double* %stddev
double* %symmat

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
--j
--Loop Bound: (0, 1024)
--Loop inc: (j + 1)
--Loop predicate: <
----array access stddev.addr j
----i
----Loop Bound: (0, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------array access data.addr ((i * 1024) + j)
------array access mean.addr j
------array access data.addr ((i * 1024) + j)
------array access mean.addr j
------array access stddev.addr j
------array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
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
------array access stddev.addr j
------array access data.addr ((i * 1024) + j)
------array access data.addr ((i * 1024) + j)
--j1
--Loop Bound: (0, 1023)
--Loop inc: (j1 + 1)
--Loop predicate: <
----array access symmat.addr ((j1 * 1024) + j1)
----j2
----Loop Bound: ((j1 + 1), 1024)
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
--array access symmat.addr 1048575

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 0 1 
Dump stride: 1 1 
init counter: 0 1 0 
Dump stride: 1 1 1 
Dump tree:
----Sample number: 1
------Sample number: 1
----Sample number: 1
------Sample number: 1
----Sample number: 1
------Sample number: 1
----Sample number: 1
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
int calAddrstddev_addr0( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrdata_addr1( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrmean_addr5( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrdata_addr2( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrmean_addr6( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrstddev_addr1( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrstddev_addr2( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrstddev_addr3( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrstddev_addr4( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrstddev_addr5( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrstddev_addr6( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrstddev_addr7( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrstddev_addr8( int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrmean_addr7( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrdata_addr3( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrdata_addr4( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrstddev_addr9( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrdata_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrdata_addr6( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr0( int j1) {
    int result = (((j1 * 1024) + j1)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr1( int j1, int j2) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrdata_addr7( int j1, int j2, int i) {
    int result = (((i * 1024) + j1)) * 8 / 64;
    return result;
}
int calAddrdata_addr8( int j1, int j2, int i) {
    int result = (((i * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr2( int j1, int j2, int i) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr3( int j1, int j2, int i) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr4( int j1, int j2) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr5( int j1, int j2) {
    int result = (((j2 * 1024) + j1)) * 8 / 64;
    return result;
}
int calAddrsymmat_addr6( ) {
    int result = (1048575) * 8 / 64;
    return result;
}
void ref_data_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_data_addr0 = -1;
    uint64_t prev_j_Start_data_addr0 = -1;
    uint64_t prev_j_End_data_addr0 = -1;
    uint64_t prev_i_Start_data_addr0 = -1;
    uint64_t prev_i_End_data_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_data_addr0 != -1) {
            if ( calAddrdata_addr0( j_Start - prev_j_Start_data_addr0 + prev_j_End_data_addr0, i_Start - prev_i_Start_data_addr0 + prev_i_End_data_addr0) == calAddrdata_addr0(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_data_addr0);
                goto EndSample;
            }
        }
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
                        prev_cnt_data_addr0 = cnt;
                        prev_j_Start_data_addr0 = j_Start;
                        prev_j_End_data_addr0 = j;
                        prev_i_Start_data_addr0 = i_Start;
                        prev_i_End_data_addr0 = i;
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
        int jLB2 = 0;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr1( j, i) == calAddrdata_addr0(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr2( j, i) == calAddrdata_addr0(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr3( i, j) == calAddrdata_addr0(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr4( i, j) == calAddrdata_addr0(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr5( i, j) == calAddrdata_addr0(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr6( i, j) == calAddrdata_addr0(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr7( j1, j2, i) == calAddrdata_addr0(j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr8( j1, j2, i) == calAddrdata_addr0(j_Start, i_Start)) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_data_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_data_addr1 = -1;
    uint64_t prev_j_Start_data_addr1 = -1;
    uint64_t prev_j_End_data_addr1 = -1;
    uint64_t prev_i_Start_data_addr1 = -1;
    uint64_t prev_i_End_data_addr1 = -1;
    uint64_t prev_cnt_data_addr2 = -1;
    uint64_t prev_j_Start_data_addr2 = -1;
    uint64_t prev_j_End_data_addr2 = -1;
    uint64_t prev_i_Start_data_addr2 = -1;
    uint64_t prev_i_End_data_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_data_addr1 != -1) {
            if ( calAddrdata_addr1( j_Start - prev_j_Start_data_addr1 + prev_j_End_data_addr1, i_Start - prev_i_Start_data_addr1 + prev_i_End_data_addr1) == calAddrdata_addr1(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_data_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr2 != -1) {
            if ( calAddrdata_addr2( j_Start - prev_j_Start_data_addr2 + prev_j_End_data_addr2, i_Start - prev_i_Start_data_addr2 + prev_i_End_data_addr2) == calAddrdata_addr1(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_data_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            if ( j == j_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr1( j, i) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr1 = cnt;
                        prev_j_Start_data_addr1 = j_Start;
                        prev_j_End_data_addr1 = j;
                        prev_i_Start_data_addr1 = i_Start;
                        prev_i_End_data_addr1 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr2( j, i) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr2 = cnt;
                        prev_j_Start_data_addr2 = j_Start;
                        prev_j_End_data_addr2 = j;
                        prev_i_Start_data_addr2 = i_Start;
                        prev_i_End_data_addr2 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr3( i, j) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr4( i, j) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr5( i, j) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr6( i, j) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr7( j1, j2, i) == calAddrdata_addr1(j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr8( j1, j2, i) == calAddrdata_addr1(j_Start, i_Start)) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_data_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_data_addr1 = -1;
    uint64_t prev_j_Start_data_addr1 = -1;
    uint64_t prev_j_End_data_addr1 = -1;
    uint64_t prev_i_Start_data_addr1 = -1;
    uint64_t prev_i_End_data_addr1 = -1;
    uint64_t prev_cnt_data_addr2 = -1;
    uint64_t prev_j_Start_data_addr2 = -1;
    uint64_t prev_j_End_data_addr2 = -1;
    uint64_t prev_i_Start_data_addr2 = -1;
    uint64_t prev_i_End_data_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_data_addr1 != -1) {
            if ( calAddrdata_addr1( j_Start - prev_j_Start_data_addr1 + prev_j_End_data_addr1, i_Start - prev_i_Start_data_addr1 + prev_i_End_data_addr1) == calAddrdata_addr2(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_data_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr2 != -1) {
            if ( calAddrdata_addr2( j_Start - prev_j_Start_data_addr2 + prev_j_End_data_addr2, i_Start - prev_i_Start_data_addr2 + prev_i_End_data_addr2) == calAddrdata_addr2(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_data_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            if ( j == j_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr1( j, i) == calAddrdata_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr1 = cnt;
                        prev_j_Start_data_addr1 = j_Start;
                        prev_j_End_data_addr1 = j;
                        prev_i_Start_data_addr1 = i_Start;
                        prev_i_End_data_addr1 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr2( j, i) == calAddrdata_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr2 = cnt;
                        prev_j_Start_data_addr2 = j_Start;
                        prev_j_End_data_addr2 = j;
                        prev_i_Start_data_addr2 = i_Start;
                        prev_i_End_data_addr2 = i;
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
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr3( i, j) == calAddrdata_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr4( i, j) == calAddrdata_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr5( i, j) == calAddrdata_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr6( i, j) == calAddrdata_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr7( j1, j2, i) == calAddrdata_addr2(j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr8( j1, j2, i) == calAddrdata_addr2(j_Start, i_Start)) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_data_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_data_addr3 = -1;
    uint64_t prev_i_Start_data_addr3 = -1;
    uint64_t prev_i_End_data_addr3 = -1;
    uint64_t prev_j_Start_data_addr3 = -1;
    uint64_t prev_j_End_data_addr3 = -1;
    uint64_t prev_cnt_data_addr4 = -1;
    uint64_t prev_i_Start_data_addr4 = -1;
    uint64_t prev_i_End_data_addr4 = -1;
    uint64_t prev_j_Start_data_addr4 = -1;
    uint64_t prev_j_End_data_addr4 = -1;
    uint64_t prev_cnt_data_addr5 = -1;
    uint64_t prev_i_Start_data_addr5 = -1;
    uint64_t prev_i_End_data_addr5 = -1;
    uint64_t prev_j_Start_data_addr5 = -1;
    uint64_t prev_j_End_data_addr5 = -1;
    uint64_t prev_cnt_data_addr6 = -1;
    uint64_t prev_i_Start_data_addr6 = -1;
    uint64_t prev_i_End_data_addr6 = -1;
    uint64_t prev_j_Start_data_addr6 = -1;
    uint64_t prev_j_End_data_addr6 = -1;
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
        if ( prev_cnt_data_addr3 != -1) {
            if ( calAddrdata_addr3( i_Start - prev_i_Start_data_addr3 + prev_i_End_data_addr3, j_Start - prev_j_Start_data_addr3 + prev_j_End_data_addr3) == calAddrdata_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr4 != -1) {
            if ( calAddrdata_addr4( i_Start - prev_i_Start_data_addr4 + prev_i_End_data_addr4, j_Start - prev_j_Start_data_addr4 + prev_j_End_data_addr4) == calAddrdata_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr5 != -1) {
            if ( calAddrdata_addr5( i_Start - prev_i_Start_data_addr5 + prev_i_End_data_addr5, j_Start - prev_j_Start_data_addr5 + prev_j_End_data_addr5) == calAddrdata_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr6 != -1) {
            if ( calAddrdata_addr6( i_Start - prev_i_Start_data_addr6 + prev_i_End_data_addr6, j_Start - prev_j_Start_data_addr6 + prev_j_End_data_addr6) == calAddrdata_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr3( i, j) == calAddrdata_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr3 = cnt;
                        prev_i_Start_data_addr3 = i_Start;
                        prev_i_End_data_addr3 = i;
                        prev_j_Start_data_addr3 = j_Start;
                        prev_j_End_data_addr3 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr4( i, j) == calAddrdata_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr4 = cnt;
                        prev_i_Start_data_addr4 = i_Start;
                        prev_i_End_data_addr4 = i;
                        prev_j_Start_data_addr4 = j_Start;
                        prev_j_End_data_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr5( i, j) == calAddrdata_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr5 = cnt;
                        prev_i_Start_data_addr5 = i_Start;
                        prev_i_End_data_addr5 = i;
                        prev_j_Start_data_addr5 = j_Start;
                        prev_j_End_data_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr6( i, j) == calAddrdata_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr6 = cnt;
                        prev_i_Start_data_addr6 = i_Start;
                        prev_i_End_data_addr6 = i;
                        prev_j_Start_data_addr6 = j_Start;
                        prev_j_End_data_addr6 = j;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr7( j1, j2, i) == calAddrdata_addr3(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr8( j1, j2, i) == calAddrdata_addr3(i_Start, j_Start)) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_data_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_data_addr3 = -1;
    uint64_t prev_i_Start_data_addr3 = -1;
    uint64_t prev_i_End_data_addr3 = -1;
    uint64_t prev_j_Start_data_addr3 = -1;
    uint64_t prev_j_End_data_addr3 = -1;
    uint64_t prev_cnt_data_addr4 = -1;
    uint64_t prev_i_Start_data_addr4 = -1;
    uint64_t prev_i_End_data_addr4 = -1;
    uint64_t prev_j_Start_data_addr4 = -1;
    uint64_t prev_j_End_data_addr4 = -1;
    uint64_t prev_cnt_data_addr5 = -1;
    uint64_t prev_i_Start_data_addr5 = -1;
    uint64_t prev_i_End_data_addr5 = -1;
    uint64_t prev_j_Start_data_addr5 = -1;
    uint64_t prev_j_End_data_addr5 = -1;
    uint64_t prev_cnt_data_addr6 = -1;
    uint64_t prev_i_Start_data_addr6 = -1;
    uint64_t prev_i_End_data_addr6 = -1;
    uint64_t prev_j_Start_data_addr6 = -1;
    uint64_t prev_j_End_data_addr6 = -1;
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
        if ( prev_cnt_data_addr3 != -1) {
            if ( calAddrdata_addr3( i_Start - prev_i_Start_data_addr3 + prev_i_End_data_addr3, j_Start - prev_j_Start_data_addr3 + prev_j_End_data_addr3) == calAddrdata_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr4 != -1) {
            if ( calAddrdata_addr4( i_Start - prev_i_Start_data_addr4 + prev_i_End_data_addr4, j_Start - prev_j_Start_data_addr4 + prev_j_End_data_addr4) == calAddrdata_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr5 != -1) {
            if ( calAddrdata_addr5( i_Start - prev_i_Start_data_addr5 + prev_i_End_data_addr5, j_Start - prev_j_Start_data_addr5 + prev_j_End_data_addr5) == calAddrdata_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr6 != -1) {
            if ( calAddrdata_addr6( i_Start - prev_i_Start_data_addr6 + prev_i_End_data_addr6, j_Start - prev_j_Start_data_addr6 + prev_j_End_data_addr6) == calAddrdata_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr3( i, j) == calAddrdata_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr3 = cnt;
                        prev_i_Start_data_addr3 = i_Start;
                        prev_i_End_data_addr3 = i;
                        prev_j_Start_data_addr3 = j_Start;
                        prev_j_End_data_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr4( i, j) == calAddrdata_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr4 = cnt;
                        prev_i_Start_data_addr4 = i_Start;
                        prev_i_End_data_addr4 = i;
                        prev_j_Start_data_addr4 = j_Start;
                        prev_j_End_data_addr4 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr5( i, j) == calAddrdata_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr5 = cnt;
                        prev_i_Start_data_addr5 = i_Start;
                        prev_i_End_data_addr5 = i;
                        prev_j_Start_data_addr5 = j_Start;
                        prev_j_End_data_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr6( i, j) == calAddrdata_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr6 = cnt;
                        prev_i_Start_data_addr6 = i_Start;
                        prev_i_End_data_addr6 = i;
                        prev_j_Start_data_addr6 = j_Start;
                        prev_j_End_data_addr6 = j;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr7( j1, j2, i) == calAddrdata_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr8( j1, j2, i) == calAddrdata_addr4(i_Start, j_Start)) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_data_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_data_addr3 = -1;
    uint64_t prev_i_Start_data_addr3 = -1;
    uint64_t prev_i_End_data_addr3 = -1;
    uint64_t prev_j_Start_data_addr3 = -1;
    uint64_t prev_j_End_data_addr3 = -1;
    uint64_t prev_cnt_data_addr4 = -1;
    uint64_t prev_i_Start_data_addr4 = -1;
    uint64_t prev_i_End_data_addr4 = -1;
    uint64_t prev_j_Start_data_addr4 = -1;
    uint64_t prev_j_End_data_addr4 = -1;
    uint64_t prev_cnt_data_addr5 = -1;
    uint64_t prev_i_Start_data_addr5 = -1;
    uint64_t prev_i_End_data_addr5 = -1;
    uint64_t prev_j_Start_data_addr5 = -1;
    uint64_t prev_j_End_data_addr5 = -1;
    uint64_t prev_cnt_data_addr6 = -1;
    uint64_t prev_i_Start_data_addr6 = -1;
    uint64_t prev_i_End_data_addr6 = -1;
    uint64_t prev_j_Start_data_addr6 = -1;
    uint64_t prev_j_End_data_addr6 = -1;
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
        if ( prev_cnt_data_addr3 != -1) {
            if ( calAddrdata_addr3( i_Start - prev_i_Start_data_addr3 + prev_i_End_data_addr3, j_Start - prev_j_Start_data_addr3 + prev_j_End_data_addr3) == calAddrdata_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr4 != -1) {
            if ( calAddrdata_addr4( i_Start - prev_i_Start_data_addr4 + prev_i_End_data_addr4, j_Start - prev_j_Start_data_addr4 + prev_j_End_data_addr4) == calAddrdata_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr5 != -1) {
            if ( calAddrdata_addr5( i_Start - prev_i_Start_data_addr5 + prev_i_End_data_addr5, j_Start - prev_j_Start_data_addr5 + prev_j_End_data_addr5) == calAddrdata_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr6 != -1) {
            if ( calAddrdata_addr6( i_Start - prev_i_Start_data_addr6 + prev_i_End_data_addr6, j_Start - prev_j_Start_data_addr6 + prev_j_End_data_addr6) == calAddrdata_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr3( i, j) == calAddrdata_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr3 = cnt;
                        prev_i_Start_data_addr3 = i_Start;
                        prev_i_End_data_addr3 = i;
                        prev_j_Start_data_addr3 = j_Start;
                        prev_j_End_data_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr4( i, j) == calAddrdata_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr4 = cnt;
                        prev_i_Start_data_addr4 = i_Start;
                        prev_i_End_data_addr4 = i;
                        prev_j_Start_data_addr4 = j_Start;
                        prev_j_End_data_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr5( i, j) == calAddrdata_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr5 = cnt;
                        prev_i_Start_data_addr5 = i_Start;
                        prev_i_End_data_addr5 = i;
                        prev_j_Start_data_addr5 = j_Start;
                        prev_j_End_data_addr5 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr6( i, j) == calAddrdata_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr6 = cnt;
                        prev_i_Start_data_addr6 = i_Start;
                        prev_i_End_data_addr6 = i;
                        prev_j_Start_data_addr6 = j_Start;
                        prev_j_End_data_addr6 = j;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr7( j1, j2, i) == calAddrdata_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr8( j1, j2, i) == calAddrdata_addr5(i_Start, j_Start)) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_data_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_data_addr3 = -1;
    uint64_t prev_i_Start_data_addr3 = -1;
    uint64_t prev_i_End_data_addr3 = -1;
    uint64_t prev_j_Start_data_addr3 = -1;
    uint64_t prev_j_End_data_addr3 = -1;
    uint64_t prev_cnt_data_addr4 = -1;
    uint64_t prev_i_Start_data_addr4 = -1;
    uint64_t prev_i_End_data_addr4 = -1;
    uint64_t prev_j_Start_data_addr4 = -1;
    uint64_t prev_j_End_data_addr4 = -1;
    uint64_t prev_cnt_data_addr5 = -1;
    uint64_t prev_i_Start_data_addr5 = -1;
    uint64_t prev_i_End_data_addr5 = -1;
    uint64_t prev_j_Start_data_addr5 = -1;
    uint64_t prev_j_End_data_addr5 = -1;
    uint64_t prev_cnt_data_addr6 = -1;
    uint64_t prev_i_Start_data_addr6 = -1;
    uint64_t prev_i_End_data_addr6 = -1;
    uint64_t prev_j_Start_data_addr6 = -1;
    uint64_t prev_j_End_data_addr6 = -1;
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
        if ( prev_cnt_data_addr3 != -1) {
            if ( calAddrdata_addr3( i_Start - prev_i_Start_data_addr3 + prev_i_End_data_addr3, j_Start - prev_j_Start_data_addr3 + prev_j_End_data_addr3) == calAddrdata_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr4 != -1) {
            if ( calAddrdata_addr4( i_Start - prev_i_Start_data_addr4 + prev_i_End_data_addr4, j_Start - prev_j_Start_data_addr4 + prev_j_End_data_addr4) == calAddrdata_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr5 != -1) {
            if ( calAddrdata_addr5( i_Start - prev_i_Start_data_addr5 + prev_i_End_data_addr5, j_Start - prev_j_Start_data_addr5 + prev_j_End_data_addr5) == calAddrdata_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr6 != -1) {
            if ( calAddrdata_addr6( i_Start - prev_i_Start_data_addr6 + prev_i_End_data_addr6, j_Start - prev_j_Start_data_addr6 + prev_j_End_data_addr6) == calAddrdata_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_data_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr3( i, j) == calAddrdata_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr3 = cnt;
                        prev_i_Start_data_addr3 = i_Start;
                        prev_i_End_data_addr3 = i;
                        prev_j_Start_data_addr3 = j_Start;
                        prev_j_End_data_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr4( i, j) == calAddrdata_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr4 = cnt;
                        prev_i_Start_data_addr4 = i_Start;
                        prev_i_End_data_addr4 = i;
                        prev_j_Start_data_addr4 = j_Start;
                        prev_j_End_data_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr5( i, j) == calAddrdata_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr5 = cnt;
                        prev_i_Start_data_addr5 = i_Start;
                        prev_i_End_data_addr5 = i;
                        prev_j_Start_data_addr5 = j_Start;
                        prev_j_End_data_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr6( i, j) == calAddrdata_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_data_addr6 = cnt;
                        prev_i_Start_data_addr6 = i_Start;
                        prev_i_End_data_addr6 = i;
                        prev_j_Start_data_addr6 = j_Start;
                        prev_j_End_data_addr6 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr7( j1, j2, i) == calAddrdata_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr8( j1, j2, i) == calAddrdata_addr6(i_Start, j_Start)) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_data_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_data_addr7 = -1;
    uint64_t prev_j1_Start_data_addr7 = -1;
    uint64_t prev_j1_End_data_addr7 = -1;
    uint64_t prev_j2_Start_data_addr7 = -1;
    uint64_t prev_j2_End_data_addr7 = -1;
    uint64_t prev_i_Start_data_addr7 = -1;
    uint64_t prev_i_End_data_addr7 = -1;
    uint64_t prev_cnt_data_addr8 = -1;
    uint64_t prev_j1_Start_data_addr8 = -1;
    uint64_t prev_j1_End_data_addr8 = -1;
    uint64_t prev_j2_Start_data_addr8 = -1;
    uint64_t prev_j2_End_data_addr8 = -1;
    uint64_t prev_i_Start_data_addr8 = -1;
    uint64_t prev_i_End_data_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_data_addr7 != -1) {
            if ( calAddrdata_addr7( j1_Start - prev_j1_Start_data_addr7 + prev_j1_End_data_addr7, j2_Start - prev_j2_Start_data_addr7 + prev_j2_End_data_addr7, i_Start - prev_i_Start_data_addr7 + prev_i_End_data_addr7) == calAddrdata_addr7(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(prev_cnt_data_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr8 != -1) {
            if ( calAddrdata_addr8( j1_Start - prev_j1_Start_data_addr8 + prev_j1_End_data_addr8, j2_Start - prev_j2_Start_data_addr8 + prev_j2_End_data_addr8, i_Start - prev_i_Start_data_addr8 + prev_i_End_data_addr8) == calAddrdata_addr7(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(prev_cnt_data_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            if ( j1 == j1_Start ) {
                j2LB7 = j2_Start;
            }
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                if ( j1 == j1_Start && j2 == j2_Start ) {
                    iLB8 = i_Start;
                }
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr7( j1, j2, i) == calAddrdata_addr7(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_data_addr7 = cnt;
                            prev_j1_Start_data_addr7 = j1_Start;
                            prev_j1_End_data_addr7 = j1;
                            prev_j2_Start_data_addr7 = j2_Start;
                            prev_j2_End_data_addr7 = j2;
                            prev_i_Start_data_addr7 = i_Start;
                            prev_i_End_data_addr7 = i;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr8( j1, j2, i) == calAddrdata_addr7(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_data_addr8 = cnt;
                            prev_j1_Start_data_addr8 = j1_Start;
                            prev_j1_End_data_addr8 = j1;
                            prev_j2_Start_data_addr8 = j2_Start;
                            prev_j2_End_data_addr8 = j2;
                            prev_i_Start_data_addr8 = i_Start;
                            prev_i_End_data_addr8 = i;
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_data_addr8() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_data_addr7 = -1;
    uint64_t prev_j1_Start_data_addr7 = -1;
    uint64_t prev_j1_End_data_addr7 = -1;
    uint64_t prev_j2_Start_data_addr7 = -1;
    uint64_t prev_j2_End_data_addr7 = -1;
    uint64_t prev_i_Start_data_addr7 = -1;
    uint64_t prev_i_End_data_addr7 = -1;
    uint64_t prev_cnt_data_addr8 = -1;
    uint64_t prev_j1_Start_data_addr8 = -1;
    uint64_t prev_j1_End_data_addr8 = -1;
    uint64_t prev_j2_Start_data_addr8 = -1;
    uint64_t prev_j2_End_data_addr8 = -1;
    uint64_t prev_i_Start_data_addr8 = -1;
    uint64_t prev_i_End_data_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_data_addr7 != -1) {
            if ( calAddrdata_addr7( j1_Start - prev_j1_Start_data_addr7 + prev_j1_End_data_addr7, j2_Start - prev_j2_Start_data_addr7 + prev_j2_End_data_addr7, i_Start - prev_i_Start_data_addr7 + prev_i_End_data_addr7) == calAddrdata_addr8(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(prev_cnt_data_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_data_addr8 != -1) {
            if ( calAddrdata_addr8( j1_Start - prev_j1_Start_data_addr8 + prev_j1_End_data_addr8, j2_Start - prev_j2_Start_data_addr8 + prev_j2_End_data_addr8, i_Start - prev_i_Start_data_addr8 + prev_i_End_data_addr8) == calAddrdata_addr8(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(prev_cnt_data_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            if ( j1 == j1_Start ) {
                j2LB7 = j2_Start;
            }
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                if ( j1 == j1_Start && j2 == j2_Start ) {
                    iLB8 = i_Start;
                }
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr7( j1, j2, i) == calAddrdata_addr8(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_data_addr7 = cnt;
                            prev_j1_Start_data_addr7 = j1_Start;
                            prev_j1_End_data_addr7 = j1;
                            prev_j2_Start_data_addr7 = j2_Start;
                            prev_j2_End_data_addr7 = j2;
                            prev_i_Start_data_addr7 = i_Start;
                            prev_i_End_data_addr7 = i;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr8( j1, j2, i) == calAddrdata_addr8(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_data_addr8 = cnt;
                            prev_j1_Start_data_addr8 = j1_Start;
                            prev_j1_End_data_addr8 = j1;
                            prev_j2_Start_data_addr8 = j2_Start;
                            prev_j2_End_data_addr8 = j2;
                            prev_i_Start_data_addr8 = i_Start;
                            prev_i_End_data_addr8 = i;
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_mean_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_mean_addr0 = -1;
    uint64_t prev_j_Start_mean_addr0 = -1;
    uint64_t prev_j_End_mean_addr0 = -1;
    uint64_t prev_cnt_mean_addr3 = -1;
    uint64_t prev_j_Start_mean_addr3 = -1;
    uint64_t prev_j_End_mean_addr3 = -1;
    uint64_t prev_cnt_mean_addr4 = -1;
    uint64_t prev_j_Start_mean_addr4 = -1;
    uint64_t prev_j_End_mean_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_mean_addr0 != -1) {
            if ( calAddrmean_addr0( j_Start - prev_j_Start_mean_addr0 + prev_j_End_mean_addr0) == calAddrmean_addr0(j_Start)) {
                rtHistoCal(prev_cnt_mean_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr3 != -1) {
            if ( calAddrmean_addr3( j_Start - prev_j_Start_mean_addr3 + prev_j_End_mean_addr3) == calAddrmean_addr0(j_Start)) {
                rtHistoCal(prev_cnt_mean_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr4 != -1) {
            if ( calAddrmean_addr4( j_Start - prev_j_Start_mean_addr4 + prev_j_End_mean_addr4) == calAddrmean_addr0(j_Start)) {
                rtHistoCal(prev_cnt_mean_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB0 = j_Start;
        for ( int j = jLB0; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr0( j) == calAddrmean_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_mean_addr0 = cnt;
                    prev_j_Start_mean_addr0 = j_Start;
                    prev_j_End_mean_addr0 = j;
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
                    prev_cnt_mean_addr3 = cnt;
                    prev_j_Start_mean_addr3 = j_Start;
                    prev_j_End_mean_addr3 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_mean_addr4 = cnt;
                    prev_j_Start_mean_addr4 = j_Start;
                    prev_j_End_mean_addr4 = j;
                    goto EndSample;
                }
            }
        }
        }
        {
        int jLB2 = 0;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( j, i) == calAddrmean_addr0(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr6( j, i) == calAddrmean_addr0(j_Start)) {
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
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr7( i, j) == calAddrmean_addr0(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_mean_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_mean_addr1 = -1;
    uint64_t prev_j_Start_mean_addr1 = -1;
    uint64_t prev_j_End_mean_addr1 = -1;
    uint64_t prev_i_Start_mean_addr1 = -1;
    uint64_t prev_i_End_mean_addr1 = -1;
    uint64_t prev_cnt_mean_addr2 = -1;
    uint64_t prev_j_Start_mean_addr2 = -1;
    uint64_t prev_j_End_mean_addr2 = -1;
    uint64_t prev_i_Start_mean_addr2 = -1;
    uint64_t prev_i_End_mean_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_mean_addr1 != -1) {
            if ( calAddrmean_addr1( j_Start - prev_j_Start_mean_addr1 + prev_j_End_mean_addr1, i_Start - prev_i_Start_mean_addr1 + prev_i_End_mean_addr1) == calAddrmean_addr1(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_mean_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr2 != -1) {
            if ( calAddrmean_addr2( j_Start - prev_j_Start_mean_addr2 + prev_j_End_mean_addr2, i_Start - prev_i_Start_mean_addr2 + prev_i_End_mean_addr2) == calAddrmean_addr1(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_mean_addr2);
                goto EndSample;
            }
        }
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
                        prev_cnt_mean_addr1 = cnt;
                        prev_j_Start_mean_addr1 = j_Start;
                        prev_j_End_mean_addr1 = j;
                        prev_i_Start_mean_addr1 = i_Start;
                        prev_i_End_mean_addr1 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_mean_addr2 = cnt;
                        prev_j_Start_mean_addr2 = j_Start;
                        prev_j_End_mean_addr2 = j;
                        prev_i_Start_mean_addr2 = i_Start;
                        prev_i_End_mean_addr2 = i;
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
        int jLB2 = 0;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( j, i) == calAddrmean_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr6( j, i) == calAddrmean_addr1(j_Start, i_Start)) {
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
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr7( i, j) == calAddrmean_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_mean_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_mean_addr1 = -1;
    uint64_t prev_j_Start_mean_addr1 = -1;
    uint64_t prev_j_End_mean_addr1 = -1;
    uint64_t prev_i_Start_mean_addr1 = -1;
    uint64_t prev_i_End_mean_addr1 = -1;
    uint64_t prev_cnt_mean_addr2 = -1;
    uint64_t prev_j_Start_mean_addr2 = -1;
    uint64_t prev_j_End_mean_addr2 = -1;
    uint64_t prev_i_Start_mean_addr2 = -1;
    uint64_t prev_i_End_mean_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_mean_addr1 != -1) {
            if ( calAddrmean_addr1( j_Start - prev_j_Start_mean_addr1 + prev_j_End_mean_addr1, i_Start - prev_i_Start_mean_addr1 + prev_i_End_mean_addr1) == calAddrmean_addr2(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_mean_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr2 != -1) {
            if ( calAddrmean_addr2( j_Start - prev_j_Start_mean_addr2 + prev_j_End_mean_addr2, i_Start - prev_i_Start_mean_addr2 + prev_i_End_mean_addr2) == calAddrmean_addr2(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_mean_addr2);
                goto EndSample;
            }
        }
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
                        prev_cnt_mean_addr1 = cnt;
                        prev_j_Start_mean_addr1 = j_Start;
                        prev_j_End_mean_addr1 = j;
                        prev_i_Start_mean_addr1 = i_Start;
                        prev_i_End_mean_addr1 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_mean_addr2 = cnt;
                        prev_j_Start_mean_addr2 = j_Start;
                        prev_j_End_mean_addr2 = j;
                        prev_i_Start_mean_addr2 = i_Start;
                        prev_i_End_mean_addr2 = i;
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
        int jLB2 = 0;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( j, i) == calAddrmean_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr6( j, i) == calAddrmean_addr2(j_Start, i_Start)) {
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
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr7( i, j) == calAddrmean_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_mean_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_mean_addr0 = -1;
    uint64_t prev_j_Start_mean_addr0 = -1;
    uint64_t prev_j_End_mean_addr0 = -1;
    uint64_t prev_cnt_mean_addr3 = -1;
    uint64_t prev_j_Start_mean_addr3 = -1;
    uint64_t prev_j_End_mean_addr3 = -1;
    uint64_t prev_cnt_mean_addr4 = -1;
    uint64_t prev_j_Start_mean_addr4 = -1;
    uint64_t prev_j_End_mean_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_mean_addr0 != -1) {
            if ( calAddrmean_addr0( j_Start - prev_j_Start_mean_addr0 + prev_j_End_mean_addr0) == calAddrmean_addr3(j_Start)) {
                rtHistoCal(prev_cnt_mean_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr3 != -1) {
            if ( calAddrmean_addr3( j_Start - prev_j_Start_mean_addr3 + prev_j_End_mean_addr3) == calAddrmean_addr3(j_Start)) {
                rtHistoCal(prev_cnt_mean_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr4 != -1) {
            if ( calAddrmean_addr4( j_Start - prev_j_Start_mean_addr4 + prev_j_End_mean_addr4) == calAddrmean_addr3(j_Start)) {
                rtHistoCal(prev_cnt_mean_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB0 = j_Start;
        for ( int j = jLB0; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr0( j) == calAddrmean_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_mean_addr0 = cnt;
                    prev_j_Start_mean_addr0 = j_Start;
                    prev_j_End_mean_addr0 = j;
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
                    prev_cnt_mean_addr3 = cnt;
                    prev_j_Start_mean_addr3 = j_Start;
                    prev_j_End_mean_addr3 = j;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_mean_addr4 = cnt;
                    prev_j_Start_mean_addr4 = j_Start;
                    prev_j_End_mean_addr4 = j;
                    goto EndSample;
                }
            }
        }
        }
        {
        int jLB2 = 0;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( j, i) == calAddrmean_addr3(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr6( j, i) == calAddrmean_addr3(j_Start)) {
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
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr7( i, j) == calAddrmean_addr3(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_mean_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_mean_addr0 = -1;
    uint64_t prev_j_Start_mean_addr0 = -1;
    uint64_t prev_j_End_mean_addr0 = -1;
    uint64_t prev_cnt_mean_addr3 = -1;
    uint64_t prev_j_Start_mean_addr3 = -1;
    uint64_t prev_j_End_mean_addr3 = -1;
    uint64_t prev_cnt_mean_addr4 = -1;
    uint64_t prev_j_Start_mean_addr4 = -1;
    uint64_t prev_j_End_mean_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_mean_addr0 != -1) {
            if ( calAddrmean_addr0( j_Start - prev_j_Start_mean_addr0 + prev_j_End_mean_addr0) == calAddrmean_addr4(j_Start)) {
                rtHistoCal(prev_cnt_mean_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr3 != -1) {
            if ( calAddrmean_addr3( j_Start - prev_j_Start_mean_addr3 + prev_j_End_mean_addr3) == calAddrmean_addr4(j_Start)) {
                rtHistoCal(prev_cnt_mean_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr4 != -1) {
            if ( calAddrmean_addr4( j_Start - prev_j_Start_mean_addr4 + prev_j_End_mean_addr4) == calAddrmean_addr4(j_Start)) {
                rtHistoCal(prev_cnt_mean_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB0 = j_Start;
        for ( int j = jLB0; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr0( j) == calAddrmean_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_mean_addr0 = cnt;
                    prev_j_Start_mean_addr0 = j_Start;
                    prev_j_End_mean_addr0 = j;
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
                    prev_cnt_mean_addr3 = cnt;
                    prev_j_Start_mean_addr3 = j_Start;
                    prev_j_End_mean_addr3 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_mean_addr4 = cnt;
                    prev_j_Start_mean_addr4 = j_Start;
                    prev_j_End_mean_addr4 = j;
                    goto EndSample;
                }
            }
            cntStart = true;
        }
        }
        {
        int jLB2 = 0;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( j, i) == calAddrmean_addr4(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr6( j, i) == calAddrmean_addr4(j_Start)) {
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
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr7( i, j) == calAddrmean_addr4(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_mean_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_mean_addr5 = -1;
    uint64_t prev_j_Start_mean_addr5 = -1;
    uint64_t prev_j_End_mean_addr5 = -1;
    uint64_t prev_i_Start_mean_addr5 = -1;
    uint64_t prev_i_End_mean_addr5 = -1;
    uint64_t prev_cnt_mean_addr6 = -1;
    uint64_t prev_j_Start_mean_addr6 = -1;
    uint64_t prev_j_End_mean_addr6 = -1;
    uint64_t prev_i_Start_mean_addr6 = -1;
    uint64_t prev_i_End_mean_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_mean_addr5 != -1) {
            if ( calAddrmean_addr5( j_Start - prev_j_Start_mean_addr5 + prev_j_End_mean_addr5, i_Start - prev_i_Start_mean_addr5 + prev_i_End_mean_addr5) == calAddrmean_addr5(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_mean_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr6 != -1) {
            if ( calAddrmean_addr6( j_Start - prev_j_Start_mean_addr6 + prev_j_End_mean_addr6, i_Start - prev_i_Start_mean_addr6 + prev_i_End_mean_addr6) == calAddrmean_addr5(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_mean_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            if ( j == j_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( j, i) == calAddrmean_addr5(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_mean_addr5 = cnt;
                        prev_j_Start_mean_addr5 = j_Start;
                        prev_j_End_mean_addr5 = j;
                        prev_i_Start_mean_addr5 = i_Start;
                        prev_i_End_mean_addr5 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr6( j, i) == calAddrmean_addr5(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_mean_addr6 = cnt;
                        prev_j_Start_mean_addr6 = j_Start;
                        prev_j_End_mean_addr6 = j;
                        prev_i_Start_mean_addr6 = i_Start;
                        prev_i_End_mean_addr6 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr7( i, j) == calAddrmean_addr5(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_mean_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_mean_addr5 = -1;
    uint64_t prev_j_Start_mean_addr5 = -1;
    uint64_t prev_j_End_mean_addr5 = -1;
    uint64_t prev_i_Start_mean_addr5 = -1;
    uint64_t prev_i_End_mean_addr5 = -1;
    uint64_t prev_cnt_mean_addr6 = -1;
    uint64_t prev_j_Start_mean_addr6 = -1;
    uint64_t prev_j_End_mean_addr6 = -1;
    uint64_t prev_i_Start_mean_addr6 = -1;
    uint64_t prev_i_End_mean_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_mean_addr5 != -1) {
            if ( calAddrmean_addr5( j_Start - prev_j_Start_mean_addr5 + prev_j_End_mean_addr5, i_Start - prev_i_Start_mean_addr5 + prev_i_End_mean_addr5) == calAddrmean_addr6(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_mean_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_mean_addr6 != -1) {
            if ( calAddrmean_addr6( j_Start - prev_j_Start_mean_addr6 + prev_j_End_mean_addr6, i_Start - prev_i_Start_mean_addr6 + prev_i_End_mean_addr6) == calAddrmean_addr6(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_mean_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            if ( j == j_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr5( j, i) == calAddrmean_addr6(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_mean_addr5 = cnt;
                        prev_j_Start_mean_addr5 = j_Start;
                        prev_j_End_mean_addr5 = j;
                        prev_i_Start_mean_addr5 = i_Start;
                        prev_i_End_mean_addr5 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr6( j, i) == calAddrmean_addr6(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_mean_addr6 = cnt;
                        prev_j_Start_mean_addr6 = j_Start;
                        prev_j_End_mean_addr6 = j;
                        prev_i_Start_mean_addr6 = i_Start;
                        prev_i_End_mean_addr6 = i;
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
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr7( i, j) == calAddrmean_addr6(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_mean_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_mean_addr7 = -1;
    uint64_t prev_i_Start_mean_addr7 = -1;
    uint64_t prev_i_End_mean_addr7 = -1;
    uint64_t prev_j_Start_mean_addr7 = -1;
    uint64_t prev_j_End_mean_addr7 = -1;
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
        if ( prev_cnt_mean_addr7 != -1) {
            if ( calAddrmean_addr7( i_Start - prev_i_Start_mean_addr7 + prev_i_End_mean_addr7, j_Start - prev_j_Start_mean_addr7 + prev_j_End_mean_addr7) == calAddrmean_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_mean_addr7);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr7( i, j) == calAddrmean_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_mean_addr7 = cnt;
                        prev_i_Start_mean_addr7 = i_Start;
                        prev_i_End_mean_addr7 = i;
                        prev_j_Start_mean_addr7 = j_Start;
                        prev_j_End_mean_addr7 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr0 = -1;
    uint64_t prev_j_Start_stddev_addr0 = -1;
    uint64_t prev_j_End_stddev_addr0 = -1;
    uint64_t prev_cnt_stddev_addr3 = -1;
    uint64_t prev_j_Start_stddev_addr3 = -1;
    uint64_t prev_j_End_stddev_addr3 = -1;
    uint64_t prev_cnt_stddev_addr4 = -1;
    uint64_t prev_j_Start_stddev_addr4 = -1;
    uint64_t prev_j_End_stddev_addr4 = -1;
    uint64_t prev_cnt_stddev_addr5 = -1;
    uint64_t prev_j_Start_stddev_addr5 = -1;
    uint64_t prev_j_End_stddev_addr5 = -1;
    uint64_t prev_cnt_stddev_addr6 = -1;
    uint64_t prev_j_Start_stddev_addr6 = -1;
    uint64_t prev_j_End_stddev_addr6 = -1;
    uint64_t prev_cnt_stddev_addr7 = -1;
    uint64_t prev_j_Start_stddev_addr7 = -1;
    uint64_t prev_j_End_stddev_addr7 = -1;
    uint64_t prev_cnt_stddev_addr8 = -1;
    uint64_t prev_j_Start_stddev_addr8 = -1;
    uint64_t prev_j_End_stddev_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_stddev_addr0 != -1) {
            if ( calAddrstddev_addr0( j_Start - prev_j_Start_stddev_addr0 + prev_j_End_stddev_addr0) == calAddrstddev_addr0(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr3 != -1) {
            if ( calAddrstddev_addr3( j_Start - prev_j_Start_stddev_addr3 + prev_j_End_stddev_addr3) == calAddrstddev_addr0(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr4 != -1) {
            if ( calAddrstddev_addr4( j_Start - prev_j_Start_stddev_addr4 + prev_j_End_stddev_addr4) == calAddrstddev_addr0(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr5 != -1) {
            if ( calAddrstddev_addr5( j_Start - prev_j_Start_stddev_addr5 + prev_j_End_stddev_addr5) == calAddrstddev_addr0(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr6 != -1) {
            if ( calAddrstddev_addr6( j_Start - prev_j_Start_stddev_addr6 + prev_j_End_stddev_addr6) == calAddrstddev_addr0(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr7 != -1) {
            if ( calAddrstddev_addr7( j_Start - prev_j_Start_stddev_addr7 + prev_j_End_stddev_addr7) == calAddrstddev_addr0(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr8 != -1) {
            if ( calAddrstddev_addr8( j_Start - prev_j_Start_stddev_addr8 + prev_j_End_stddev_addr8) == calAddrstddev_addr0(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr0( j) == calAddrstddev_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr0 = cnt;
                    prev_j_Start_stddev_addr0 = j_Start;
                    prev_j_End_stddev_addr0 = j;
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr1( j, i) == calAddrstddev_addr0(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr2( j, i) == calAddrstddev_addr0(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr3( j) == calAddrstddev_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr3 = cnt;
                    prev_j_Start_stddev_addr3 = j_Start;
                    prev_j_End_stddev_addr3 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr4( j) == calAddrstddev_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr4 = cnt;
                    prev_j_Start_stddev_addr4 = j_Start;
                    prev_j_End_stddev_addr4 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr5( j) == calAddrstddev_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr5 = cnt;
                    prev_j_Start_stddev_addr5 = j_Start;
                    prev_j_End_stddev_addr5 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr6 = cnt;
                    prev_j_Start_stddev_addr6 = j_Start;
                    prev_j_End_stddev_addr6 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr7( j) == calAddrstddev_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr7 = cnt;
                    prev_j_Start_stddev_addr7 = j_Start;
                    prev_j_End_stddev_addr7 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr8( j) == calAddrstddev_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr8 = cnt;
                    prev_j_Start_stddev_addr8 = j_Start;
                    prev_j_End_stddev_addr8 = j;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr0(j_Start)) {
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr1 = -1;
    uint64_t prev_j_Start_stddev_addr1 = -1;
    uint64_t prev_j_End_stddev_addr1 = -1;
    uint64_t prev_i_Start_stddev_addr1 = -1;
    uint64_t prev_i_End_stddev_addr1 = -1;
    uint64_t prev_cnt_stddev_addr2 = -1;
    uint64_t prev_j_Start_stddev_addr2 = -1;
    uint64_t prev_j_End_stddev_addr2 = -1;
    uint64_t prev_i_Start_stddev_addr2 = -1;
    uint64_t prev_i_End_stddev_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_stddev_addr1 != -1) {
            if ( calAddrstddev_addr1( j_Start - prev_j_Start_stddev_addr1 + prev_j_End_stddev_addr1, i_Start - prev_i_Start_stddev_addr1 + prev_i_End_stddev_addr1) == calAddrstddev_addr1(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_stddev_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr2 != -1) {
            if ( calAddrstddev_addr2( j_Start - prev_j_Start_stddev_addr2 + prev_j_End_stddev_addr2, i_Start - prev_i_Start_stddev_addr2 + prev_i_End_stddev_addr2) == calAddrstddev_addr1(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_stddev_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr0( j) == calAddrstddev_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB3 = 0;
            if ( j == j_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr1( j, i) == calAddrstddev_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_stddev_addr1 = cnt;
                        prev_j_Start_stddev_addr1 = j_Start;
                        prev_j_End_stddev_addr1 = j;
                        prev_i_Start_stddev_addr1 = i_Start;
                        prev_i_End_stddev_addr1 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr2( j, i) == calAddrstddev_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_stddev_addr2 = cnt;
                        prev_j_Start_stddev_addr2 = j_Start;
                        prev_j_End_stddev_addr2 = j;
                        prev_i_Start_stddev_addr2 = i_Start;
                        prev_i_End_stddev_addr2 = i;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr3( j) == calAddrstddev_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr4( j) == calAddrstddev_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr5( j) == calAddrstddev_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr7( j) == calAddrstddev_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr8( j) == calAddrstddev_addr1(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr1(j_Start, i_Start)) {
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr1 = -1;
    uint64_t prev_j_Start_stddev_addr1 = -1;
    uint64_t prev_j_End_stddev_addr1 = -1;
    uint64_t prev_i_Start_stddev_addr1 = -1;
    uint64_t prev_i_End_stddev_addr1 = -1;
    uint64_t prev_cnt_stddev_addr2 = -1;
    uint64_t prev_j_Start_stddev_addr2 = -1;
    uint64_t prev_j_End_stddev_addr2 = -1;
    uint64_t prev_i_Start_stddev_addr2 = -1;
    uint64_t prev_i_End_stddev_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_stddev_addr1 != -1) {
            if ( calAddrstddev_addr1( j_Start - prev_j_Start_stddev_addr1 + prev_j_End_stddev_addr1, i_Start - prev_i_Start_stddev_addr1 + prev_i_End_stddev_addr1) == calAddrstddev_addr2(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_stddev_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr2 != -1) {
            if ( calAddrstddev_addr2( j_Start - prev_j_Start_stddev_addr2 + prev_j_End_stddev_addr2, i_Start - prev_i_Start_stddev_addr2 + prev_i_End_stddev_addr2) == calAddrstddev_addr2(j_Start, i_Start)) {
                rtHistoCal(prev_cnt_stddev_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr0( j) == calAddrstddev_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB3 = 0;
            if ( j == j_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr1( j, i) == calAddrstddev_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_stddev_addr1 = cnt;
                        prev_j_Start_stddev_addr1 = j_Start;
                        prev_j_End_stddev_addr1 = j;
                        prev_i_Start_stddev_addr1 = i_Start;
                        prev_i_End_stddev_addr1 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr2( j, i) == calAddrstddev_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_stddev_addr2 = cnt;
                        prev_j_Start_stddev_addr2 = j_Start;
                        prev_j_End_stddev_addr2 = j;
                        prev_i_Start_stddev_addr2 = i_Start;
                        prev_i_End_stddev_addr2 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr3( j) == calAddrstddev_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr4( j) == calAddrstddev_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr5( j) == calAddrstddev_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr7( j) == calAddrstddev_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr8( j) == calAddrstddev_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr2(j_Start, i_Start)) {
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr0 = -1;
    uint64_t prev_j_Start_stddev_addr0 = -1;
    uint64_t prev_j_End_stddev_addr0 = -1;
    uint64_t prev_cnt_stddev_addr3 = -1;
    uint64_t prev_j_Start_stddev_addr3 = -1;
    uint64_t prev_j_End_stddev_addr3 = -1;
    uint64_t prev_cnt_stddev_addr4 = -1;
    uint64_t prev_j_Start_stddev_addr4 = -1;
    uint64_t prev_j_End_stddev_addr4 = -1;
    uint64_t prev_cnt_stddev_addr5 = -1;
    uint64_t prev_j_Start_stddev_addr5 = -1;
    uint64_t prev_j_End_stddev_addr5 = -1;
    uint64_t prev_cnt_stddev_addr6 = -1;
    uint64_t prev_j_Start_stddev_addr6 = -1;
    uint64_t prev_j_End_stddev_addr6 = -1;
    uint64_t prev_cnt_stddev_addr7 = -1;
    uint64_t prev_j_Start_stddev_addr7 = -1;
    uint64_t prev_j_End_stddev_addr7 = -1;
    uint64_t prev_cnt_stddev_addr8 = -1;
    uint64_t prev_j_Start_stddev_addr8 = -1;
    uint64_t prev_j_End_stddev_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_stddev_addr0 != -1) {
            if ( calAddrstddev_addr0( j_Start - prev_j_Start_stddev_addr0 + prev_j_End_stddev_addr0) == calAddrstddev_addr3(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr3 != -1) {
            if ( calAddrstddev_addr3( j_Start - prev_j_Start_stddev_addr3 + prev_j_End_stddev_addr3) == calAddrstddev_addr3(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr4 != -1) {
            if ( calAddrstddev_addr4( j_Start - prev_j_Start_stddev_addr4 + prev_j_End_stddev_addr4) == calAddrstddev_addr3(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr5 != -1) {
            if ( calAddrstddev_addr5( j_Start - prev_j_Start_stddev_addr5 + prev_j_End_stddev_addr5) == calAddrstddev_addr3(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr6 != -1) {
            if ( calAddrstddev_addr6( j_Start - prev_j_Start_stddev_addr6 + prev_j_End_stddev_addr6) == calAddrstddev_addr3(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr7 != -1) {
            if ( calAddrstddev_addr7( j_Start - prev_j_Start_stddev_addr7 + prev_j_End_stddev_addr7) == calAddrstddev_addr3(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr8 != -1) {
            if ( calAddrstddev_addr8( j_Start - prev_j_Start_stddev_addr8 + prev_j_End_stddev_addr8) == calAddrstddev_addr3(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr0( j) == calAddrstddev_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr0 = cnt;
                    prev_j_Start_stddev_addr0 = j_Start;
                    prev_j_End_stddev_addr0 = j;
                    goto EndSample;
                }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr1( j, i) == calAddrstddev_addr3(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr2( j, i) == calAddrstddev_addr3(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr3( j) == calAddrstddev_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr3 = cnt;
                    prev_j_Start_stddev_addr3 = j_Start;
                    prev_j_End_stddev_addr3 = j;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr4( j) == calAddrstddev_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr4 = cnt;
                    prev_j_Start_stddev_addr4 = j_Start;
                    prev_j_End_stddev_addr4 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr5( j) == calAddrstddev_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr5 = cnt;
                    prev_j_Start_stddev_addr5 = j_Start;
                    prev_j_End_stddev_addr5 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr6 = cnt;
                    prev_j_Start_stddev_addr6 = j_Start;
                    prev_j_End_stddev_addr6 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr7( j) == calAddrstddev_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr7 = cnt;
                    prev_j_Start_stddev_addr7 = j_Start;
                    prev_j_End_stddev_addr7 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr8( j) == calAddrstddev_addr3(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr8 = cnt;
                    prev_j_Start_stddev_addr8 = j_Start;
                    prev_j_End_stddev_addr8 = j;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr3(j_Start)) {
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr0 = -1;
    uint64_t prev_j_Start_stddev_addr0 = -1;
    uint64_t prev_j_End_stddev_addr0 = -1;
    uint64_t prev_cnt_stddev_addr3 = -1;
    uint64_t prev_j_Start_stddev_addr3 = -1;
    uint64_t prev_j_End_stddev_addr3 = -1;
    uint64_t prev_cnt_stddev_addr4 = -1;
    uint64_t prev_j_Start_stddev_addr4 = -1;
    uint64_t prev_j_End_stddev_addr4 = -1;
    uint64_t prev_cnt_stddev_addr5 = -1;
    uint64_t prev_j_Start_stddev_addr5 = -1;
    uint64_t prev_j_End_stddev_addr5 = -1;
    uint64_t prev_cnt_stddev_addr6 = -1;
    uint64_t prev_j_Start_stddev_addr6 = -1;
    uint64_t prev_j_End_stddev_addr6 = -1;
    uint64_t prev_cnt_stddev_addr7 = -1;
    uint64_t prev_j_Start_stddev_addr7 = -1;
    uint64_t prev_j_End_stddev_addr7 = -1;
    uint64_t prev_cnt_stddev_addr8 = -1;
    uint64_t prev_j_Start_stddev_addr8 = -1;
    uint64_t prev_j_End_stddev_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_stddev_addr0 != -1) {
            if ( calAddrstddev_addr0( j_Start - prev_j_Start_stddev_addr0 + prev_j_End_stddev_addr0) == calAddrstddev_addr4(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr3 != -1) {
            if ( calAddrstddev_addr3( j_Start - prev_j_Start_stddev_addr3 + prev_j_End_stddev_addr3) == calAddrstddev_addr4(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr4 != -1) {
            if ( calAddrstddev_addr4( j_Start - prev_j_Start_stddev_addr4 + prev_j_End_stddev_addr4) == calAddrstddev_addr4(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr5 != -1) {
            if ( calAddrstddev_addr5( j_Start - prev_j_Start_stddev_addr5 + prev_j_End_stddev_addr5) == calAddrstddev_addr4(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr6 != -1) {
            if ( calAddrstddev_addr6( j_Start - prev_j_Start_stddev_addr6 + prev_j_End_stddev_addr6) == calAddrstddev_addr4(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr7 != -1) {
            if ( calAddrstddev_addr7( j_Start - prev_j_Start_stddev_addr7 + prev_j_End_stddev_addr7) == calAddrstddev_addr4(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr8 != -1) {
            if ( calAddrstddev_addr8( j_Start - prev_j_Start_stddev_addr8 + prev_j_End_stddev_addr8) == calAddrstddev_addr4(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr0( j) == calAddrstddev_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr0 = cnt;
                    prev_j_Start_stddev_addr0 = j_Start;
                    prev_j_End_stddev_addr0 = j;
                    goto EndSample;
                }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr1( j, i) == calAddrstddev_addr4(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr2( j, i) == calAddrstddev_addr4(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr3( j) == calAddrstddev_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr3 = cnt;
                    prev_j_Start_stddev_addr3 = j_Start;
                    prev_j_End_stddev_addr3 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr4( j) == calAddrstddev_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr4 = cnt;
                    prev_j_Start_stddev_addr4 = j_Start;
                    prev_j_End_stddev_addr4 = j;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr5( j) == calAddrstddev_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr5 = cnt;
                    prev_j_Start_stddev_addr5 = j_Start;
                    prev_j_End_stddev_addr5 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr6 = cnt;
                    prev_j_Start_stddev_addr6 = j_Start;
                    prev_j_End_stddev_addr6 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr7( j) == calAddrstddev_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr7 = cnt;
                    prev_j_Start_stddev_addr7 = j_Start;
                    prev_j_End_stddev_addr7 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr8( j) == calAddrstddev_addr4(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr8 = cnt;
                    prev_j_Start_stddev_addr8 = j_Start;
                    prev_j_End_stddev_addr8 = j;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr4(j_Start)) {
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr0 = -1;
    uint64_t prev_j_Start_stddev_addr0 = -1;
    uint64_t prev_j_End_stddev_addr0 = -1;
    uint64_t prev_cnt_stddev_addr3 = -1;
    uint64_t prev_j_Start_stddev_addr3 = -1;
    uint64_t prev_j_End_stddev_addr3 = -1;
    uint64_t prev_cnt_stddev_addr4 = -1;
    uint64_t prev_j_Start_stddev_addr4 = -1;
    uint64_t prev_j_End_stddev_addr4 = -1;
    uint64_t prev_cnt_stddev_addr5 = -1;
    uint64_t prev_j_Start_stddev_addr5 = -1;
    uint64_t prev_j_End_stddev_addr5 = -1;
    uint64_t prev_cnt_stddev_addr6 = -1;
    uint64_t prev_j_Start_stddev_addr6 = -1;
    uint64_t prev_j_End_stddev_addr6 = -1;
    uint64_t prev_cnt_stddev_addr7 = -1;
    uint64_t prev_j_Start_stddev_addr7 = -1;
    uint64_t prev_j_End_stddev_addr7 = -1;
    uint64_t prev_cnt_stddev_addr8 = -1;
    uint64_t prev_j_Start_stddev_addr8 = -1;
    uint64_t prev_j_End_stddev_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_stddev_addr0 != -1) {
            if ( calAddrstddev_addr0( j_Start - prev_j_Start_stddev_addr0 + prev_j_End_stddev_addr0) == calAddrstddev_addr5(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr3 != -1) {
            if ( calAddrstddev_addr3( j_Start - prev_j_Start_stddev_addr3 + prev_j_End_stddev_addr3) == calAddrstddev_addr5(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr4 != -1) {
            if ( calAddrstddev_addr4( j_Start - prev_j_Start_stddev_addr4 + prev_j_End_stddev_addr4) == calAddrstddev_addr5(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr5 != -1) {
            if ( calAddrstddev_addr5( j_Start - prev_j_Start_stddev_addr5 + prev_j_End_stddev_addr5) == calAddrstddev_addr5(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr6 != -1) {
            if ( calAddrstddev_addr6( j_Start - prev_j_Start_stddev_addr6 + prev_j_End_stddev_addr6) == calAddrstddev_addr5(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr7 != -1) {
            if ( calAddrstddev_addr7( j_Start - prev_j_Start_stddev_addr7 + prev_j_End_stddev_addr7) == calAddrstddev_addr5(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr8 != -1) {
            if ( calAddrstddev_addr8( j_Start - prev_j_Start_stddev_addr8 + prev_j_End_stddev_addr8) == calAddrstddev_addr5(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr0( j) == calAddrstddev_addr5(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr0 = cnt;
                    prev_j_Start_stddev_addr0 = j_Start;
                    prev_j_End_stddev_addr0 = j;
                    goto EndSample;
                }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr1( j, i) == calAddrstddev_addr5(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr2( j, i) == calAddrstddev_addr5(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr3( j) == calAddrstddev_addr5(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr3 = cnt;
                    prev_j_Start_stddev_addr3 = j_Start;
                    prev_j_End_stddev_addr3 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr4( j) == calAddrstddev_addr5(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr4 = cnt;
                    prev_j_Start_stddev_addr4 = j_Start;
                    prev_j_End_stddev_addr4 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr5( j) == calAddrstddev_addr5(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr5 = cnt;
                    prev_j_Start_stddev_addr5 = j_Start;
                    prev_j_End_stddev_addr5 = j;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr5(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr6 = cnt;
                    prev_j_Start_stddev_addr6 = j_Start;
                    prev_j_End_stddev_addr6 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr7( j) == calAddrstddev_addr5(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr7 = cnt;
                    prev_j_Start_stddev_addr7 = j_Start;
                    prev_j_End_stddev_addr7 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr8( j) == calAddrstddev_addr5(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr8 = cnt;
                    prev_j_Start_stddev_addr8 = j_Start;
                    prev_j_End_stddev_addr8 = j;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr5(j_Start)) {
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr0 = -1;
    uint64_t prev_j_Start_stddev_addr0 = -1;
    uint64_t prev_j_End_stddev_addr0 = -1;
    uint64_t prev_cnt_stddev_addr3 = -1;
    uint64_t prev_j_Start_stddev_addr3 = -1;
    uint64_t prev_j_End_stddev_addr3 = -1;
    uint64_t prev_cnt_stddev_addr4 = -1;
    uint64_t prev_j_Start_stddev_addr4 = -1;
    uint64_t prev_j_End_stddev_addr4 = -1;
    uint64_t prev_cnt_stddev_addr5 = -1;
    uint64_t prev_j_Start_stddev_addr5 = -1;
    uint64_t prev_j_End_stddev_addr5 = -1;
    uint64_t prev_cnt_stddev_addr6 = -1;
    uint64_t prev_j_Start_stddev_addr6 = -1;
    uint64_t prev_j_End_stddev_addr6 = -1;
    uint64_t prev_cnt_stddev_addr7 = -1;
    uint64_t prev_j_Start_stddev_addr7 = -1;
    uint64_t prev_j_End_stddev_addr7 = -1;
    uint64_t prev_cnt_stddev_addr8 = -1;
    uint64_t prev_j_Start_stddev_addr8 = -1;
    uint64_t prev_j_End_stddev_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_stddev_addr0 != -1) {
            if ( calAddrstddev_addr0( j_Start - prev_j_Start_stddev_addr0 + prev_j_End_stddev_addr0) == calAddrstddev_addr6(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr3 != -1) {
            if ( calAddrstddev_addr3( j_Start - prev_j_Start_stddev_addr3 + prev_j_End_stddev_addr3) == calAddrstddev_addr6(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr4 != -1) {
            if ( calAddrstddev_addr4( j_Start - prev_j_Start_stddev_addr4 + prev_j_End_stddev_addr4) == calAddrstddev_addr6(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr5 != -1) {
            if ( calAddrstddev_addr5( j_Start - prev_j_Start_stddev_addr5 + prev_j_End_stddev_addr5) == calAddrstddev_addr6(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr6 != -1) {
            if ( calAddrstddev_addr6( j_Start - prev_j_Start_stddev_addr6 + prev_j_End_stddev_addr6) == calAddrstddev_addr6(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr7 != -1) {
            if ( calAddrstddev_addr7( j_Start - prev_j_Start_stddev_addr7 + prev_j_End_stddev_addr7) == calAddrstddev_addr6(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr8 != -1) {
            if ( calAddrstddev_addr8( j_Start - prev_j_Start_stddev_addr8 + prev_j_End_stddev_addr8) == calAddrstddev_addr6(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr0( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr0 = cnt;
                    prev_j_Start_stddev_addr0 = j_Start;
                    prev_j_End_stddev_addr0 = j;
                    goto EndSample;
                }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr1( j, i) == calAddrstddev_addr6(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr2( j, i) == calAddrstddev_addr6(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr3( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr3 = cnt;
                    prev_j_Start_stddev_addr3 = j_Start;
                    prev_j_End_stddev_addr3 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr4( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr4 = cnt;
                    prev_j_Start_stddev_addr4 = j_Start;
                    prev_j_End_stddev_addr4 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr5( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr5 = cnt;
                    prev_j_Start_stddev_addr5 = j_Start;
                    prev_j_End_stddev_addr5 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr6 = cnt;
                    prev_j_Start_stddev_addr6 = j_Start;
                    prev_j_End_stddev_addr6 = j;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr7( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr7 = cnt;
                    prev_j_Start_stddev_addr7 = j_Start;
                    prev_j_End_stddev_addr7 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr8( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr8 = cnt;
                    prev_j_Start_stddev_addr8 = j_Start;
                    prev_j_End_stddev_addr8 = j;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr6(j_Start)) {
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr0 = -1;
    uint64_t prev_j_Start_stddev_addr0 = -1;
    uint64_t prev_j_End_stddev_addr0 = -1;
    uint64_t prev_cnt_stddev_addr3 = -1;
    uint64_t prev_j_Start_stddev_addr3 = -1;
    uint64_t prev_j_End_stddev_addr3 = -1;
    uint64_t prev_cnt_stddev_addr4 = -1;
    uint64_t prev_j_Start_stddev_addr4 = -1;
    uint64_t prev_j_End_stddev_addr4 = -1;
    uint64_t prev_cnt_stddev_addr5 = -1;
    uint64_t prev_j_Start_stddev_addr5 = -1;
    uint64_t prev_j_End_stddev_addr5 = -1;
    uint64_t prev_cnt_stddev_addr6 = -1;
    uint64_t prev_j_Start_stddev_addr6 = -1;
    uint64_t prev_j_End_stddev_addr6 = -1;
    uint64_t prev_cnt_stddev_addr7 = -1;
    uint64_t prev_j_Start_stddev_addr7 = -1;
    uint64_t prev_j_End_stddev_addr7 = -1;
    uint64_t prev_cnt_stddev_addr8 = -1;
    uint64_t prev_j_Start_stddev_addr8 = -1;
    uint64_t prev_j_End_stddev_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_stddev_addr0 != -1) {
            if ( calAddrstddev_addr0( j_Start - prev_j_Start_stddev_addr0 + prev_j_End_stddev_addr0) == calAddrstddev_addr7(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr3 != -1) {
            if ( calAddrstddev_addr3( j_Start - prev_j_Start_stddev_addr3 + prev_j_End_stddev_addr3) == calAddrstddev_addr7(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr4 != -1) {
            if ( calAddrstddev_addr4( j_Start - prev_j_Start_stddev_addr4 + prev_j_End_stddev_addr4) == calAddrstddev_addr7(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr5 != -1) {
            if ( calAddrstddev_addr5( j_Start - prev_j_Start_stddev_addr5 + prev_j_End_stddev_addr5) == calAddrstddev_addr7(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr6 != -1) {
            if ( calAddrstddev_addr6( j_Start - prev_j_Start_stddev_addr6 + prev_j_End_stddev_addr6) == calAddrstddev_addr7(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr7 != -1) {
            if ( calAddrstddev_addr7( j_Start - prev_j_Start_stddev_addr7 + prev_j_End_stddev_addr7) == calAddrstddev_addr7(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr8 != -1) {
            if ( calAddrstddev_addr8( j_Start - prev_j_Start_stddev_addr8 + prev_j_End_stddev_addr8) == calAddrstddev_addr7(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr0( j) == calAddrstddev_addr7(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr0 = cnt;
                    prev_j_Start_stddev_addr0 = j_Start;
                    prev_j_End_stddev_addr0 = j;
                    goto EndSample;
                }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr1( j, i) == calAddrstddev_addr7(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr2( j, i) == calAddrstddev_addr7(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr3( j) == calAddrstddev_addr7(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr3 = cnt;
                    prev_j_Start_stddev_addr3 = j_Start;
                    prev_j_End_stddev_addr3 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr4( j) == calAddrstddev_addr7(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr4 = cnt;
                    prev_j_Start_stddev_addr4 = j_Start;
                    prev_j_End_stddev_addr4 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr5( j) == calAddrstddev_addr7(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr5 = cnt;
                    prev_j_Start_stddev_addr5 = j_Start;
                    prev_j_End_stddev_addr5 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr7(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr6 = cnt;
                    prev_j_Start_stddev_addr6 = j_Start;
                    prev_j_End_stddev_addr6 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr7( j) == calAddrstddev_addr7(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr7 = cnt;
                    prev_j_Start_stddev_addr7 = j_Start;
                    prev_j_End_stddev_addr7 = j;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr8( j) == calAddrstddev_addr7(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr8 = cnt;
                    prev_j_Start_stddev_addr8 = j_Start;
                    prev_j_End_stddev_addr8 = j;
                    goto EndSample;
                }
            }
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr7(j_Start)) {
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr8() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr0 = -1;
    uint64_t prev_j_Start_stddev_addr0 = -1;
    uint64_t prev_j_End_stddev_addr0 = -1;
    uint64_t prev_cnt_stddev_addr3 = -1;
    uint64_t prev_j_Start_stddev_addr3 = -1;
    uint64_t prev_j_End_stddev_addr3 = -1;
    uint64_t prev_cnt_stddev_addr4 = -1;
    uint64_t prev_j_Start_stddev_addr4 = -1;
    uint64_t prev_j_End_stddev_addr4 = -1;
    uint64_t prev_cnt_stddev_addr5 = -1;
    uint64_t prev_j_Start_stddev_addr5 = -1;
    uint64_t prev_j_End_stddev_addr5 = -1;
    uint64_t prev_cnt_stddev_addr6 = -1;
    uint64_t prev_j_Start_stddev_addr6 = -1;
    uint64_t prev_j_End_stddev_addr6 = -1;
    uint64_t prev_cnt_stddev_addr7 = -1;
    uint64_t prev_j_Start_stddev_addr7 = -1;
    uint64_t prev_j_End_stddev_addr7 = -1;
    uint64_t prev_cnt_stddev_addr8 = -1;
    uint64_t prev_j_Start_stddev_addr8 = -1;
    uint64_t prev_j_End_stddev_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_stddev_addr0 != -1) {
            if ( calAddrstddev_addr0( j_Start - prev_j_Start_stddev_addr0 + prev_j_End_stddev_addr0) == calAddrstddev_addr8(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr3 != -1) {
            if ( calAddrstddev_addr3( j_Start - prev_j_Start_stddev_addr3 + prev_j_End_stddev_addr3) == calAddrstddev_addr8(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr4 != -1) {
            if ( calAddrstddev_addr4( j_Start - prev_j_Start_stddev_addr4 + prev_j_End_stddev_addr4) == calAddrstddev_addr8(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr5 != -1) {
            if ( calAddrstddev_addr5( j_Start - prev_j_Start_stddev_addr5 + prev_j_End_stddev_addr5) == calAddrstddev_addr8(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr6 != -1) {
            if ( calAddrstddev_addr6( j_Start - prev_j_Start_stddev_addr6 + prev_j_End_stddev_addr6) == calAddrstddev_addr8(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr7 != -1) {
            if ( calAddrstddev_addr7( j_Start - prev_j_Start_stddev_addr7 + prev_j_End_stddev_addr7) == calAddrstddev_addr8(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_stddev_addr8 != -1) {
            if ( calAddrstddev_addr8( j_Start - prev_j_Start_stddev_addr8 + prev_j_End_stddev_addr8) == calAddrstddev_addr8(j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr0( j) == calAddrstddev_addr8(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr0 = cnt;
                    prev_j_Start_stddev_addr0 = j_Start;
                    prev_j_End_stddev_addr0 = j;
                    goto EndSample;
                }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr1( j, i) == calAddrstddev_addr8(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr2( j, i) == calAddrstddev_addr8(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr3( j) == calAddrstddev_addr8(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr3 = cnt;
                    prev_j_Start_stddev_addr3 = j_Start;
                    prev_j_End_stddev_addr3 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr4( j) == calAddrstddev_addr8(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr4 = cnt;
                    prev_j_Start_stddev_addr4 = j_Start;
                    prev_j_End_stddev_addr4 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr5( j) == calAddrstddev_addr8(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr5 = cnt;
                    prev_j_Start_stddev_addr5 = j_Start;
                    prev_j_End_stddev_addr5 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr8(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr6 = cnt;
                    prev_j_Start_stddev_addr6 = j_Start;
                    prev_j_End_stddev_addr6 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr7( j) == calAddrstddev_addr8(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr7 = cnt;
                    prev_j_Start_stddev_addr7 = j_Start;
                    prev_j_End_stddev_addr7 = j;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr8( j) == calAddrstddev_addr8(j_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_stddev_addr8 = cnt;
                    prev_j_Start_stddev_addr8 = j_Start;
                    prev_j_End_stddev_addr8 = j;
                    goto EndSample;
                }
            }
            cntStart = true;
        }
        }
        {
        int iLB4 = 0;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr8(j_Start)) {
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_stddev_addr9() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_stddev_addr9 = -1;
    uint64_t prev_i_Start_stddev_addr9 = -1;
    uint64_t prev_i_End_stddev_addr9 = -1;
    uint64_t prev_j_Start_stddev_addr9 = -1;
    uint64_t prev_j_End_stddev_addr9 = -1;
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
        if ( prev_cnt_stddev_addr9 != -1) {
            if ( calAddrstddev_addr9( i_Start - prev_i_Start_stddev_addr9 + prev_i_End_stddev_addr9, j_Start - prev_j_Start_stddev_addr9 + prev_j_End_stddev_addr9) == calAddrstddev_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_stddev_addr9);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr9( i, j) == calAddrstddev_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_stddev_addr9 = cnt;
                        prev_i_Start_stddev_addr9 = i_Start;
                        prev_i_End_stddev_addr9 = i;
                        prev_j_Start_stddev_addr9 = j_Start;
                        prev_j_End_stddev_addr9 = j;
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
        int j1LB6 = 0;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) cnt++;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) cnt++;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_symmat_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_symmat_addr0 = -1;
    uint64_t prev_j1_Start_symmat_addr0 = -1;
    uint64_t prev_j1_End_symmat_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_symmat_addr0 != -1) {
            if ( calAddrsymmat_addr0( j1_Start - prev_j1_Start_symmat_addr0 + prev_j1_End_symmat_addr0) == calAddrsymmat_addr0(j1_Start)) {
                rtHistoCal(prev_cnt_symmat_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr0( j1) == calAddrsymmat_addr0(j1_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_symmat_addr0 = cnt;
                    prev_j1_Start_symmat_addr0 = j1_Start;
                    prev_j1_End_symmat_addr0 = j1;
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr1( j1, j2) == calAddrsymmat_addr0(j1_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr0(j1_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr3( j1, j2, i) == calAddrsymmat_addr0(j1_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr0(j1_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr5( j1, j2) == calAddrsymmat_addr0(j1_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        if (cntStart == true) {
            cnt++;
            if ( calAddrsymmat_addr6( ) == calAddrsymmat_addr0(j1_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_symmat_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_symmat_addr1 = -1;
    uint64_t prev_j1_Start_symmat_addr1 = -1;
    uint64_t prev_j1_End_symmat_addr1 = -1;
    uint64_t prev_j2_Start_symmat_addr1 = -1;
    uint64_t prev_j2_End_symmat_addr1 = -1;
    uint64_t prev_cnt_symmat_addr4 = -1;
    uint64_t prev_j1_Start_symmat_addr4 = -1;
    uint64_t prev_j1_End_symmat_addr4 = -1;
    uint64_t prev_j2_Start_symmat_addr4 = -1;
    uint64_t prev_j2_End_symmat_addr4 = -1;
    uint64_t prev_cnt_symmat_addr5 = -1;
    uint64_t prev_j1_Start_symmat_addr5 = -1;
    uint64_t prev_j1_End_symmat_addr5 = -1;
    uint64_t prev_j2_Start_symmat_addr5 = -1;
    uint64_t prev_j2_End_symmat_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_symmat_addr1 != -1) {
            if ( calAddrsymmat_addr1( j1_Start - prev_j1_Start_symmat_addr1 + prev_j1_End_symmat_addr1, j2_Start - prev_j2_Start_symmat_addr1 + prev_j2_End_symmat_addr1) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                rtHistoCal(prev_cnt_symmat_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_symmat_addr4 != -1) {
            if ( calAddrsymmat_addr4( j1_Start - prev_j1_Start_symmat_addr4 + prev_j1_End_symmat_addr4, j2_Start - prev_j2_Start_symmat_addr4 + prev_j2_End_symmat_addr4) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                rtHistoCal(prev_cnt_symmat_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_symmat_addr5 != -1) {
            if ( calAddrsymmat_addr5( j1_Start - prev_j1_Start_symmat_addr5 + prev_j1_End_symmat_addr5, j2_Start - prev_j2_Start_symmat_addr5 + prev_j2_End_symmat_addr5) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                rtHistoCal(prev_cnt_symmat_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr0( j1) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int j2LB7 = (j1 + 1);
            if ( j1 == j1_Start ) {
                j2LB7 = j2_Start;
            }
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr1( j1, j2) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_symmat_addr1 = cnt;
                        prev_j1_Start_symmat_addr1 = j1_Start;
                        prev_j1_End_symmat_addr1 = j1;
                        prev_j2_Start_symmat_addr1 = j2_Start;
                        prev_j2_End_symmat_addr1 = j2;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr3( j1, j2, i) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_symmat_addr4 = cnt;
                        prev_j1_Start_symmat_addr4 = j1_Start;
                        prev_j1_End_symmat_addr4 = j1;
                        prev_j2_Start_symmat_addr4 = j2_Start;
                        prev_j2_End_symmat_addr4 = j2;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr5( j1, j2) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_symmat_addr5 = cnt;
                        prev_j1_Start_symmat_addr5 = j1_Start;
                        prev_j1_End_symmat_addr5 = j1;
                        prev_j2_Start_symmat_addr5 = j2_Start;
                        prev_j2_End_symmat_addr5 = j2;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        if (cntStart == true) {
            cnt++;
            if ( calAddrsymmat_addr6( ) == calAddrsymmat_addr1(j1_Start, j2_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_symmat_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_symmat_addr2 = -1;
    uint64_t prev_j1_Start_symmat_addr2 = -1;
    uint64_t prev_j1_End_symmat_addr2 = -1;
    uint64_t prev_j2_Start_symmat_addr2 = -1;
    uint64_t prev_j2_End_symmat_addr2 = -1;
    uint64_t prev_i_Start_symmat_addr2 = -1;
    uint64_t prev_i_End_symmat_addr2 = -1;
    uint64_t prev_cnt_symmat_addr3 = -1;
    uint64_t prev_j1_Start_symmat_addr3 = -1;
    uint64_t prev_j1_End_symmat_addr3 = -1;
    uint64_t prev_j2_Start_symmat_addr3 = -1;
    uint64_t prev_j2_End_symmat_addr3 = -1;
    uint64_t prev_i_Start_symmat_addr3 = -1;
    uint64_t prev_i_End_symmat_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_symmat_addr2 != -1) {
            if ( calAddrsymmat_addr2( j1_Start - prev_j1_Start_symmat_addr2 + prev_j1_End_symmat_addr2, j2_Start - prev_j2_Start_symmat_addr2 + prev_j2_End_symmat_addr2, i_Start - prev_i_Start_symmat_addr2 + prev_i_End_symmat_addr2) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(prev_cnt_symmat_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_symmat_addr3 != -1) {
            if ( calAddrsymmat_addr3( j1_Start - prev_j1_Start_symmat_addr3 + prev_j1_End_symmat_addr3, j2_Start - prev_j2_Start_symmat_addr3 + prev_j2_End_symmat_addr3, i_Start - prev_i_Start_symmat_addr3 + prev_i_End_symmat_addr3) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(prev_cnt_symmat_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr0( j1) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int j2LB7 = (j1 + 1);
            if ( j1 == j1_Start ) {
                j2LB7 = j2_Start;
            }
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr1( j1, j2) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB8 = 0;
                if ( j1 == j1_Start && j2 == j2_Start ) {
                    iLB8 = i_Start;
                }
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_symmat_addr2 = cnt;
                            prev_j1_Start_symmat_addr2 = j1_Start;
                            prev_j1_End_symmat_addr2 = j1;
                            prev_j2_Start_symmat_addr2 = j2_Start;
                            prev_j2_End_symmat_addr2 = j2;
                            prev_i_Start_symmat_addr2 = i_Start;
                            prev_i_End_symmat_addr2 = i;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr3( j1, j2, i) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_symmat_addr3 = cnt;
                            prev_j1_Start_symmat_addr3 = j1_Start;
                            prev_j1_End_symmat_addr3 = j1;
                            prev_j2_Start_symmat_addr3 = j2_Start;
                            prev_j2_End_symmat_addr3 = j2;
                            prev_i_Start_symmat_addr3 = i_Start;
                            prev_i_End_symmat_addr3 = i;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr5( j1, j2) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        if (cntStart == true) {
            cnt++;
            if ( calAddrsymmat_addr6( ) == calAddrsymmat_addr2(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_symmat_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_symmat_addr2 = -1;
    uint64_t prev_j1_Start_symmat_addr2 = -1;
    uint64_t prev_j1_End_symmat_addr2 = -1;
    uint64_t prev_j2_Start_symmat_addr2 = -1;
    uint64_t prev_j2_End_symmat_addr2 = -1;
    uint64_t prev_i_Start_symmat_addr2 = -1;
    uint64_t prev_i_End_symmat_addr2 = -1;
    uint64_t prev_cnt_symmat_addr3 = -1;
    uint64_t prev_j1_Start_symmat_addr3 = -1;
    uint64_t prev_j1_End_symmat_addr3 = -1;
    uint64_t prev_j2_Start_symmat_addr3 = -1;
    uint64_t prev_j2_End_symmat_addr3 = -1;
    uint64_t prev_i_Start_symmat_addr3 = -1;
    uint64_t prev_i_End_symmat_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_symmat_addr2 != -1) {
            if ( calAddrsymmat_addr2( j1_Start - prev_j1_Start_symmat_addr2 + prev_j1_End_symmat_addr2, j2_Start - prev_j2_Start_symmat_addr2 + prev_j2_End_symmat_addr2, i_Start - prev_i_Start_symmat_addr2 + prev_i_End_symmat_addr2) == calAddrsymmat_addr3(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(prev_cnt_symmat_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_symmat_addr3 != -1) {
            if ( calAddrsymmat_addr3( j1_Start - prev_j1_Start_symmat_addr3 + prev_j1_End_symmat_addr3, j2_Start - prev_j2_Start_symmat_addr3 + prev_j2_End_symmat_addr3, i_Start - prev_i_Start_symmat_addr3 + prev_i_End_symmat_addr3) == calAddrsymmat_addr3(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(prev_cnt_symmat_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr0( j1) == calAddrsymmat_addr3(j1_Start, j2_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int j2LB7 = (j1 + 1);
            if ( j1 == j1_Start ) {
                j2LB7 = j2_Start;
            }
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr1( j1, j2) == calAddrsymmat_addr3(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int iLB8 = 0;
                if ( j1 == j1_Start && j2 == j2_Start ) {
                    iLB8 = i_Start;
                }
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr3(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_symmat_addr2 = cnt;
                            prev_j1_Start_symmat_addr2 = j1_Start;
                            prev_j1_End_symmat_addr2 = j1;
                            prev_j2_Start_symmat_addr2 = j2_Start;
                            prev_j2_End_symmat_addr2 = j2;
                            prev_i_Start_symmat_addr2 = i_Start;
                            prev_i_End_symmat_addr2 = i;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr3( j1, j2, i) == calAddrsymmat_addr3(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_symmat_addr3 = cnt;
                            prev_j1_Start_symmat_addr3 = j1_Start;
                            prev_j1_End_symmat_addr3 = j1;
                            prev_j2_Start_symmat_addr3 = j2_Start;
                            prev_j2_End_symmat_addr3 = j2;
                            prev_i_Start_symmat_addr3 = i_Start;
                            prev_i_End_symmat_addr3 = i;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr3(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr5( j1, j2) == calAddrsymmat_addr3(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        if (cntStart == true) {
            cnt++;
            if ( calAddrsymmat_addr6( ) == calAddrsymmat_addr3(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_symmat_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_symmat_addr1 = -1;
    uint64_t prev_j1_Start_symmat_addr1 = -1;
    uint64_t prev_j1_End_symmat_addr1 = -1;
    uint64_t prev_j2_Start_symmat_addr1 = -1;
    uint64_t prev_j2_End_symmat_addr1 = -1;
    uint64_t prev_cnt_symmat_addr4 = -1;
    uint64_t prev_j1_Start_symmat_addr4 = -1;
    uint64_t prev_j1_End_symmat_addr4 = -1;
    uint64_t prev_j2_Start_symmat_addr4 = -1;
    uint64_t prev_j2_End_symmat_addr4 = -1;
    uint64_t prev_cnt_symmat_addr5 = -1;
    uint64_t prev_j1_Start_symmat_addr5 = -1;
    uint64_t prev_j1_End_symmat_addr5 = -1;
    uint64_t prev_j2_Start_symmat_addr5 = -1;
    uint64_t prev_j2_End_symmat_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_symmat_addr1 != -1) {
            if ( calAddrsymmat_addr1( j1_Start - prev_j1_Start_symmat_addr1 + prev_j1_End_symmat_addr1, j2_Start - prev_j2_Start_symmat_addr1 + prev_j2_End_symmat_addr1) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                rtHistoCal(prev_cnt_symmat_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_symmat_addr4 != -1) {
            if ( calAddrsymmat_addr4( j1_Start - prev_j1_Start_symmat_addr4 + prev_j1_End_symmat_addr4, j2_Start - prev_j2_Start_symmat_addr4 + prev_j2_End_symmat_addr4) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                rtHistoCal(prev_cnt_symmat_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_symmat_addr5 != -1) {
            if ( calAddrsymmat_addr5( j1_Start - prev_j1_Start_symmat_addr5 + prev_j1_End_symmat_addr5, j2_Start - prev_j2_Start_symmat_addr5 + prev_j2_End_symmat_addr5) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                rtHistoCal(prev_cnt_symmat_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr0( j1) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int j2LB7 = (j1 + 1);
            if ( j1 == j1_Start ) {
                j2LB7 = j2_Start;
            }
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr1( j1, j2) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_symmat_addr1 = cnt;
                        prev_j1_Start_symmat_addr1 = j1_Start;
                        prev_j1_End_symmat_addr1 = j1;
                        prev_j2_Start_symmat_addr1 = j2_Start;
                        prev_j2_End_symmat_addr1 = j2;
                        goto EndSample;
                    }
                }
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr3( j1, j2, i) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_symmat_addr4 = cnt;
                        prev_j1_Start_symmat_addr4 = j1_Start;
                        prev_j1_End_symmat_addr4 = j1;
                        prev_j2_Start_symmat_addr4 = j2_Start;
                        prev_j2_End_symmat_addr4 = j2;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr5( j1, j2) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_symmat_addr5 = cnt;
                        prev_j1_Start_symmat_addr5 = j1_Start;
                        prev_j1_End_symmat_addr5 = j1;
                        prev_j2_Start_symmat_addr5 = j2_Start;
                        prev_j2_End_symmat_addr5 = j2;
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        if (cntStart == true) {
            cnt++;
            if ( calAddrsymmat_addr6( ) == calAddrsymmat_addr4(j1_Start, j2_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_symmat_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_symmat_addr1 = -1;
    uint64_t prev_j1_Start_symmat_addr1 = -1;
    uint64_t prev_j1_End_symmat_addr1 = -1;
    uint64_t prev_j2_Start_symmat_addr1 = -1;
    uint64_t prev_j2_End_symmat_addr1 = -1;
    uint64_t prev_cnt_symmat_addr4 = -1;
    uint64_t prev_j1_Start_symmat_addr4 = -1;
    uint64_t prev_j1_End_symmat_addr4 = -1;
    uint64_t prev_j2_Start_symmat_addr4 = -1;
    uint64_t prev_j2_End_symmat_addr4 = -1;
    uint64_t prev_cnt_symmat_addr5 = -1;
    uint64_t prev_j1_Start_symmat_addr5 = -1;
    uint64_t prev_j1_End_symmat_addr5 = -1;
    uint64_t prev_j2_Start_symmat_addr5 = -1;
    uint64_t prev_j2_End_symmat_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 0;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_symmat_addr1 != -1) {
            if ( calAddrsymmat_addr1( j1_Start - prev_j1_Start_symmat_addr1 + prev_j1_End_symmat_addr1, j2_Start - prev_j2_Start_symmat_addr1 + prev_j2_End_symmat_addr1) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                rtHistoCal(prev_cnt_symmat_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_symmat_addr4 != -1) {
            if ( calAddrsymmat_addr4( j1_Start - prev_j1_Start_symmat_addr4 + prev_j1_End_symmat_addr4, j2_Start - prev_j2_Start_symmat_addr4 + prev_j2_End_symmat_addr4) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                rtHistoCal(prev_cnt_symmat_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_symmat_addr5 != -1) {
            if ( calAddrsymmat_addr5( j1_Start - prev_j1_Start_symmat_addr5 + prev_j1_End_symmat_addr5, j2_Start - prev_j2_Start_symmat_addr5 + prev_j2_End_symmat_addr5) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                rtHistoCal(prev_cnt_symmat_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr0( j1) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int j2LB7 = (j1 + 1);
            if ( j1 == j1_Start ) {
                j2LB7 = j2_Start;
            }
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr1( j1, j2) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_symmat_addr1 = cnt;
                        prev_j1_Start_symmat_addr1 = j1_Start;
                        prev_j1_End_symmat_addr1 = j1;
                        prev_j2_Start_symmat_addr1 = j2_Start;
                        prev_j2_End_symmat_addr1 = j2;
                        goto EndSample;
                    }
                }
                {
                int iLB8 = 0;
                for ( int i = iLB8; i < 1024; i++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr2( j1, j2, i) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr3( j1, j2, i) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr4( j1, j2) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_symmat_addr4 = cnt;
                        prev_j1_Start_symmat_addr4 = j1_Start;
                        prev_j1_End_symmat_addr4 = j1;
                        prev_j2_Start_symmat_addr4 = j2_Start;
                        prev_j2_End_symmat_addr4 = j2;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr5( j1, j2) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_symmat_addr5 = cnt;
                        prev_j1_Start_symmat_addr5 = j1_Start;
                        prev_j1_End_symmat_addr5 = j1;
                        prev_j2_Start_symmat_addr5 = j2_Start;
                        prev_j2_End_symmat_addr5 = j2;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        if (cntStart == true) {
            cnt++;
            if ( calAddrsymmat_addr6( ) == calAddrsymmat_addr5(j1_Start, j2_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_symmat_addr6() {
    /* Generating search reuse init code (same loop) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
    /* Generating reuse search code */

    s++;
    }
}
int main() {
    ref_data_addr0();
    ref_data_addr1();
    ref_data_addr2();
    ref_data_addr3();
    ref_data_addr4();
    ref_data_addr5();
    ref_data_addr6();
    ref_data_addr7();
    ref_data_addr8();
    ref_mean_addr0();
    ref_mean_addr1();
    ref_mean_addr2();
    ref_mean_addr3();
    ref_mean_addr4();
    ref_mean_addr5();
    ref_mean_addr6();
    ref_mean_addr7();
    ref_stddev_addr0();
    ref_stddev_addr1();
    ref_stddev_addr2();
    ref_stddev_addr3();
    ref_stddev_addr4();
    ref_stddev_addr5();
    ref_stddev_addr6();
    ref_stddev_addr7();
    ref_stddev_addr8();
    ref_stddev_addr9();
    ref_symmat_addr0();
    ref_symmat_addr1();
    ref_symmat_addr2();
    ref_symmat_addr3();
    ref_symmat_addr4();
    ref_symmat_addr5();
    ref_symmat_addr6();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: correlation */ 
