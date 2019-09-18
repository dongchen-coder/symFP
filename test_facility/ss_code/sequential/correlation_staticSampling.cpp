
 /* Start to analysis array index
Array index info: Total number of references: 34
stddev.addr j
mean.addr j
data.addr ((i * 1024) + j)
mean.addr j
stddev.addr j
stddev.addr j
mean.addr j
mean.addr j
mean.addr j
mean.addr j
stddev.addr j
stddev.addr j
data.addr ((i * 1024) + j)
mean.addr j
data.addr ((i * 1024) + j)
symmat.addr ((j1 * 1024) + j2)
data.addr ((i * 1024) + j1)
data.addr ((i * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
stddev.addr j
stddev.addr j
symmat.addr ((j1 * 1024) + j2)
stddev.addr j
stddev.addr j
mean.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
stddev.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
symmat.addr ((j1 * 1024) + j1)
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
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 52
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
/* mean_addr j 0 */
int calAddrmean_addr0( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* data_addr ((i * 1024) + j) 1 */
int calAddrdata_addr1( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* mean_addr j 2 */
int calAddrmean_addr2( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
/* mean_addr j 3 */
int calAddrmean_addr3( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
/* mean_addr j 4 */
int calAddrmean_addr4( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* mean_addr j 5 */
int calAddrmean_addr5( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* stddev_addr j 6 */
int calAddrstddev_addr6( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* data_addr ((i * 1024) + j) 7 */
int calAddrdata_addr7( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* mean_addr j 8 */
int calAddrmean_addr8( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
/* data_addr ((i * 1024) + j) 9 */
int calAddrdata_addr9( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* mean_addr j 10 */
int calAddrmean_addr10( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
/* stddev_addr j 11 */
int calAddrstddev_addr11( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
/* stddev_addr j 12 */
int calAddrstddev_addr12( int j, int i) {
    int result = (j) * 8 / 64;
    return result;
}
/* stddev_addr j 13 */
int calAddrstddev_addr13( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* stddev_addr j 14 */
int calAddrstddev_addr14( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* stddev_addr j 15 */
int calAddrstddev_addr15( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* stddev_addr j 16 */
int calAddrstddev_addr16( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* stddev_addr j 17 */
int calAddrstddev_addr17( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* stddev_addr j 18 */
int calAddrstddev_addr18( int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* mean_addr j 19 */
int calAddrmean_addr19( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* data_addr ((i * 1024) + j) 20 */
int calAddrdata_addr20( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* data_addr ((i * 1024) + j) 21 */
int calAddrdata_addr21( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* stddev_addr j 22 */
int calAddrstddev_addr22( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* data_addr ((i * 1024) + j) 23 */
int calAddrdata_addr23( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* data_addr ((i * 1024) + j) 24 */
int calAddrdata_addr24( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* symmat_addr ((j1 * 1024) + j1) 25 */
int calAddrsymmat_addr25( int j1) {
    int result = (((j1 * 1024) + j1)) * 8 / 64;
    return result;
}
/* symmat_addr ((j1 * 1024) + j2) 26 */
int calAddrsymmat_addr26( int j1, int j2) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
/* data_addr ((i * 1024) + j1) 27 */
int calAddrdata_addr27( int j1, int j2, int i) {
    int result = (((i * 1024) + j1)) * 8 / 64;
    return result;
}
/* data_addr ((i * 1024) + j2) 28 */
int calAddrdata_addr28( int j1, int j2, int i) {
    int result = (((i * 1024) + j2)) * 8 / 64;
    return result;
}
/* symmat_addr ((j1 * 1024) + j2) 29 */
int calAddrsymmat_addr29( int j1, int j2, int i) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
/* symmat_addr ((j1 * 1024) + j2) 30 */
int calAddrsymmat_addr30( int j1, int j2, int i) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
/* symmat_addr ((j1 * 1024) + j2) 31 */
int calAddrsymmat_addr31( int j1, int j2) {
    int result = (((j1 * 1024) + j2)) * 8 / 64;
    return result;
}
/* symmat_addr ((j2 * 1024) + j1) 32 */
int calAddrsymmat_addr32( int j1, int j2) {
    int result = (((j2 * 1024) + j1)) * 8 / 64;
    return result;
}
/* symmat_addr 1048575 33 */
int calAddrsymmat_addr33( ) {
    int result = (1048575) * 8 / 64;
    return result;
}
void ref_stddev_addr6() {
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
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
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
                    if ( calAddrstddev_addr11( j, i) == calAddrstddev_addr6(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr12( j, i) == calAddrstddev_addr6(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr13( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr14( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr15( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr16( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr17( j) == calAddrstddev_addr6(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr18( j) == calAddrstddev_addr6(j_Start)) {
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr6(j_Start)) {
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
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr0(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr3( j, i) == calAddrmean_addr0(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr0(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr5( j) == calAddrmean_addr0(j_Start)) {
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
                    if ( calAddrmean_addr8( j, i) == calAddrmean_addr0(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr10( j, i) == calAddrmean_addr0(j_Start)) {
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
                    if ( calAddrmean_addr19( i, j) == calAddrmean_addr0(j_Start)) {
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
void ref_data_addr1() {
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
                    if ( calAddrdata_addr1( j, i) == calAddrdata_addr1(j_Start, i_Start)) {
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
        int jLB2 = 0;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) cnt++;
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr7( j, i) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr9( j, i) == calAddrdata_addr1(j_Start, i_Start)) {
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
                    if ( calAddrdata_addr20( i, j) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr21( i, j) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr23( i, j) == calAddrdata_addr1(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr24( i, j) == calAddrdata_addr1(j_Start, i_Start)) {
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
                        if ( calAddrdata_addr27( j1, j2, i) == calAddrdata_addr1(j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr28( j1, j2, i) == calAddrdata_addr1(j_Start, i_Start)) {
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
void ref_mean_addr10() {
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
                    if ( calAddrmean_addr8( j, i) == calAddrmean_addr10(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr10( j, i) == calAddrmean_addr10(j_Start, i_Start)) {
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
                    if ( calAddrmean_addr19( i, j) == calAddrmean_addr10(j_Start, i_Start)) {
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
void ref_stddev_addr11() {
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
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr11(j_Start, i_Start)) {
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
                    if ( calAddrstddev_addr11( j, i) == calAddrstddev_addr11(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr12( j, i) == calAddrstddev_addr11(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr13( j) == calAddrstddev_addr11(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr14( j) == calAddrstddev_addr11(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr15( j) == calAddrstddev_addr11(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr16( j) == calAddrstddev_addr11(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr17( j) == calAddrstddev_addr11(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr18( j) == calAddrstddev_addr11(j_Start, i_Start)) {
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr11(j_Start, i_Start)) {
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
void ref_stddev_addr12() {
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
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr12(j_Start, i_Start)) {
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
                    if ( calAddrstddev_addr11( j, i) == calAddrstddev_addr12(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr12( j, i) == calAddrstddev_addr12(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr13( j) == calAddrstddev_addr12(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr14( j) == calAddrstddev_addr12(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr15( j) == calAddrstddev_addr12(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr16( j) == calAddrstddev_addr12(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr17( j) == calAddrstddev_addr12(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr18( j) == calAddrstddev_addr12(j_Start, i_Start)) {
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr12(j_Start, i_Start)) {
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
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr3( j, i) == calAddrmean_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr2(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr5( j) == calAddrmean_addr2(j_Start, i_Start)) {
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
                    if ( calAddrmean_addr8( j, i) == calAddrmean_addr2(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr10( j, i) == calAddrmean_addr2(j_Start, i_Start)) {
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
                    if ( calAddrmean_addr19( i, j) == calAddrmean_addr2(j_Start, i_Start)) {
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
                if ( calAddrmean_addr0( j) == calAddrmean_addr3(j_Start, i_Start)) {
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
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr3(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr3( j, i) == calAddrmean_addr3(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr3(j_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr5( j) == calAddrmean_addr3(j_Start, i_Start)) {
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
                    if ( calAddrmean_addr8( j, i) == calAddrmean_addr3(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr10( j, i) == calAddrmean_addr3(j_Start, i_Start)) {
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
                    if ( calAddrmean_addr19( i, j) == calAddrmean_addr3(j_Start, i_Start)) {
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
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr4(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr3( j, i) == calAddrmean_addr4(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
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
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr5( j) == calAddrmean_addr4(j_Start)) {
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
                    if ( calAddrmean_addr8( j, i) == calAddrmean_addr4(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr10( j, i) == calAddrmean_addr4(j_Start)) {
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
                    if ( calAddrmean_addr19( i, j) == calAddrmean_addr4(j_Start)) {
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
                if ( calAddrmean_addr0( j) == calAddrmean_addr5(j_Start)) {
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
                    if ( calAddrmean_addr2( j, i) == calAddrmean_addr5(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr3( j, i) == calAddrmean_addr5(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr4( j) == calAddrmean_addr5(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrmean_addr5( j) == calAddrmean_addr5(j_Start)) {
                    rtHistoCal(cnt);
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
                    if ( calAddrmean_addr8( j, i) == calAddrmean_addr5(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr10( j, i) == calAddrmean_addr5(j_Start)) {
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
                    if ( calAddrmean_addr19( i, j) == calAddrmean_addr5(j_Start)) {
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
void ref_stddev_addr17() {
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
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr17(j_Start)) {
                    rtHistoCal(cnt);
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
                    if ( calAddrstddev_addr11( j, i) == calAddrstddev_addr17(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr12( j, i) == calAddrstddev_addr17(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr13( j) == calAddrstddev_addr17(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr14( j) == calAddrstddev_addr17(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr15( j) == calAddrstddev_addr17(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr16( j) == calAddrstddev_addr17(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr17( j) == calAddrstddev_addr17(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr18( j) == calAddrstddev_addr17(j_Start)) {
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr17(j_Start)) {
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
void ref_stddev_addr18() {
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
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr18(j_Start)) {
                    rtHistoCal(cnt);
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
                    if ( calAddrstddev_addr11( j, i) == calAddrstddev_addr18(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr12( j, i) == calAddrstddev_addr18(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr13( j) == calAddrstddev_addr18(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr14( j) == calAddrstddev_addr18(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr15( j) == calAddrstddev_addr18(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr16( j) == calAddrstddev_addr18(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr17( j) == calAddrstddev_addr18(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr18( j) == calAddrstddev_addr18(j_Start)) {
                    rtHistoCal(cnt);
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr18(j_Start)) {
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
void ref_data_addr7() {
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
                    if ( calAddrdata_addr7( j, i) == calAddrdata_addr7(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr9( j, i) == calAddrdata_addr7(j_Start, i_Start)) {
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
                    if ( calAddrdata_addr20( i, j) == calAddrdata_addr7(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr21( i, j) == calAddrdata_addr7(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr23( i, j) == calAddrdata_addr7(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr24( i, j) == calAddrdata_addr7(j_Start, i_Start)) {
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
                        if ( calAddrdata_addr27( j1, j2, i) == calAddrdata_addr7(j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr28( j1, j2, i) == calAddrdata_addr7(j_Start, i_Start)) {
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
void ref_mean_addr8() {
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
                    if ( calAddrmean_addr8( j, i) == calAddrmean_addr8(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrmean_addr10( j, i) == calAddrmean_addr8(j_Start, i_Start)) {
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
                    if ( calAddrmean_addr19( i, j) == calAddrmean_addr8(j_Start, i_Start)) {
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
void ref_data_addr9() {
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
                    if ( calAddrdata_addr7( j, i) == calAddrdata_addr9(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr9( j, i) == calAddrdata_addr9(j_Start, i_Start)) {
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
                    if ( calAddrdata_addr20( i, j) == calAddrdata_addr9(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr21( i, j) == calAddrdata_addr9(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr23( i, j) == calAddrdata_addr9(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr24( i, j) == calAddrdata_addr9(j_Start, i_Start)) {
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
                        if ( calAddrdata_addr27( j1, j2, i) == calAddrdata_addr9(j_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr28( j1, j2, i) == calAddrdata_addr9(j_Start, i_Start)) {
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
void ref_symmat_addr26() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr25( j1) == calAddrsymmat_addr26(j1_Start, j2_Start)) {
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
                    if ( calAddrsymmat_addr26( j1, j2) == calAddrsymmat_addr26(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
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
                        if ( calAddrsymmat_addr29( j1, j2, i) == calAddrsymmat_addr26(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr30( j1, j2, i) == calAddrsymmat_addr26(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr31( j1, j2) == calAddrsymmat_addr26(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr32( j1, j2) == calAddrsymmat_addr26(j1_Start, j2_Start)) {
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
            if ( calAddrsymmat_addr33( ) == calAddrsymmat_addr26(j1_Start, j2_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_data_addr27() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
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
                        if ( calAddrdata_addr27( j1, j2, i) == calAddrdata_addr27(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr28( j1, j2, i) == calAddrdata_addr27(j1_Start, j2_Start, i_Start)) {
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
void ref_data_addr28() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
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
                        if ( calAddrdata_addr27( j1, j2, i) == calAddrdata_addr28(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr28( j1, j2, i) == calAddrdata_addr28(j1_Start, j2_Start, i_Start)) {
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
        if (cntStart == true) cnt++;
EndSample:
        s++;
        }
}
void ref_symmat_addr29() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
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

        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr25( j1) == calAddrsymmat_addr29(j1_Start, j2_Start, i_Start)) {
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
                    if ( calAddrsymmat_addr26( j1, j2) == calAddrsymmat_addr29(j1_Start, j2_Start, i_Start)) {
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
                        if ( calAddrsymmat_addr29( j1, j2, i) == calAddrsymmat_addr29(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr30( j1, j2, i) == calAddrsymmat_addr29(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr31( j1, j2) == calAddrsymmat_addr29(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr32( j1, j2) == calAddrsymmat_addr29(j1_Start, j2_Start, i_Start)) {
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
            if ( calAddrsymmat_addr33( ) == calAddrsymmat_addr29(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_symmat_addr30() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 536;) {
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

        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr25( j1) == calAddrsymmat_addr30(j1_Start, j2_Start, i_Start)) {
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
                    if ( calAddrsymmat_addr26( j1, j2) == calAddrsymmat_addr30(j1_Start, j2_Start, i_Start)) {
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
                        if ( calAddrsymmat_addr29( j1, j2, i) == calAddrsymmat_addr30(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr30( j1, j2, i) == calAddrsymmat_addr30(j1_Start, j2_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr31( j1, j2) == calAddrsymmat_addr30(j1_Start, j2_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr32( j1, j2) == calAddrsymmat_addr30(j1_Start, j2_Start, i_Start)) {
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
            if ( calAddrsymmat_addr33( ) == calAddrsymmat_addr30(j1_Start, j2_Start, i_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_stddev_addr13() {
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
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr13(j_Start)) {
                    rtHistoCal(cnt);
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
                    if ( calAddrstddev_addr11( j, i) == calAddrstddev_addr13(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr12( j, i) == calAddrstddev_addr13(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr13( j) == calAddrstddev_addr13(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr14( j) == calAddrstddev_addr13(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr15( j) == calAddrstddev_addr13(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr16( j) == calAddrstddev_addr13(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr17( j) == calAddrstddev_addr13(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr18( j) == calAddrstddev_addr13(j_Start)) {
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr13(j_Start)) {
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
void ref_stddev_addr14() {
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
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr14(j_Start)) {
                    rtHistoCal(cnt);
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
                    if ( calAddrstddev_addr11( j, i) == calAddrstddev_addr14(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr12( j, i) == calAddrstddev_addr14(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr13( j) == calAddrstddev_addr14(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr14( j) == calAddrstddev_addr14(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr15( j) == calAddrstddev_addr14(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr16( j) == calAddrstddev_addr14(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr17( j) == calAddrstddev_addr14(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr18( j) == calAddrstddev_addr14(j_Start)) {
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr14(j_Start)) {
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
void ref_symmat_addr31() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr25( j1) == calAddrsymmat_addr31(j1_Start, j2_Start)) {
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
                    if ( calAddrsymmat_addr26( j1, j2) == calAddrsymmat_addr31(j1_Start, j2_Start)) {
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
                        if ( calAddrsymmat_addr29( j1, j2, i) == calAddrsymmat_addr31(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr30( j1, j2, i) == calAddrsymmat_addr31(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr31( j1, j2) == calAddrsymmat_addr31(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr32( j1, j2) == calAddrsymmat_addr31(j1_Start, j2_Start)) {
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
            if ( calAddrsymmat_addr33( ) == calAddrsymmat_addr31(j1_Start, j2_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_stddev_addr15() {
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
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr15(j_Start)) {
                    rtHistoCal(cnt);
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
                    if ( calAddrstddev_addr11( j, i) == calAddrstddev_addr15(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr12( j, i) == calAddrstddev_addr15(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr13( j) == calAddrstddev_addr15(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr14( j) == calAddrstddev_addr15(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr15( j) == calAddrstddev_addr15(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr16( j) == calAddrstddev_addr15(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr17( j) == calAddrstddev_addr15(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr18( j) == calAddrstddev_addr15(j_Start)) {
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr15(j_Start)) {
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
void ref_stddev_addr16() {
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
        int jLB2 = j_Start;
        for ( int j = jLB2; j < 1024; j++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr6( j) == calAddrstddev_addr16(j_Start)) {
                    rtHistoCal(cnt);
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
                    if ( calAddrstddev_addr11( j, i) == calAddrstddev_addr16(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrstddev_addr12( j, i) == calAddrstddev_addr16(j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr13( j) == calAddrstddev_addr16(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr14( j) == calAddrstddev_addr16(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr15( j) == calAddrstddev_addr16(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr16( j) == calAddrstddev_addr16(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr17( j) == calAddrstddev_addr16(j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrstddev_addr18( j) == calAddrstddev_addr16(j_Start)) {
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr16(j_Start)) {
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
void ref_mean_addr19() {
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
                    if ( calAddrmean_addr19( i, j) == calAddrmean_addr19(i_Start, j_Start)) {
                        rtHistoCal(cnt);
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
void ref_data_addr20() {
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
                    if ( calAddrdata_addr20( i, j) == calAddrdata_addr20(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr21( i, j) == calAddrdata_addr20(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr23( i, j) == calAddrdata_addr20(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr24( i, j) == calAddrdata_addr20(i_Start, j_Start)) {
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
                        if ( calAddrdata_addr27( j1, j2, i) == calAddrdata_addr20(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr28( j1, j2, i) == calAddrdata_addr20(i_Start, j_Start)) {
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
void ref_data_addr21() {
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
                    if ( calAddrdata_addr20( i, j) == calAddrdata_addr21(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr21( i, j) == calAddrdata_addr21(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr23( i, j) == calAddrdata_addr21(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr24( i, j) == calAddrdata_addr21(i_Start, j_Start)) {
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
                        if ( calAddrdata_addr27( j1, j2, i) == calAddrdata_addr21(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr28( j1, j2, i) == calAddrdata_addr21(i_Start, j_Start)) {
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
void ref_stddev_addr22() {
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
                    if ( calAddrstddev_addr22( i, j) == calAddrstddev_addr22(i_Start, j_Start)) {
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
void ref_data_addr23() {
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
                    if ( calAddrdata_addr20( i, j) == calAddrdata_addr23(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr21( i, j) == calAddrdata_addr23(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr23( i, j) == calAddrdata_addr23(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr24( i, j) == calAddrdata_addr23(i_Start, j_Start)) {
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
                        if ( calAddrdata_addr27( j1, j2, i) == calAddrdata_addr23(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr28( j1, j2, i) == calAddrdata_addr23(i_Start, j_Start)) {
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
void ref_data_addr24() {
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
                    if ( calAddrdata_addr20( i, j) == calAddrdata_addr24(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr21( i, j) == calAddrdata_addr24(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr23( i, j) == calAddrdata_addr24(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrdata_addr24( i, j) == calAddrdata_addr24(i_Start, j_Start)) {
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
                        if ( calAddrdata_addr27( j1, j2, i) == calAddrdata_addr24(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrdata_addr28( j1, j2, i) == calAddrdata_addr24(i_Start, j_Start)) {
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
void ref_symmat_addr25() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(j1_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr25( j1) == calAddrsymmat_addr25(j1_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int j2LB7 = (j1 + 1);
            for ( int j2 = j2LB7; j2 < 1024; j2++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr26( j1, j2) == calAddrsymmat_addr25(j1_Start)) {
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
                        if ( calAddrsymmat_addr29( j1, j2, i) == calAddrsymmat_addr25(j1_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr30( j1, j2, i) == calAddrsymmat_addr25(j1_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr31( j1, j2) == calAddrsymmat_addr25(j1_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr32( j1, j2) == calAddrsymmat_addr25(j1_Start)) {
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
            if ( calAddrsymmat_addr33( ) == calAddrsymmat_addr25(j1_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_symmat_addr32() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int j1_Start = rand() % (1023 - 0) + 0;
        if ( (1024 - (j1_Start + 1)) == 0) goto SAMPLE;
        int j2_Start = rand() % (1024 - (j1_Start + 1)) + (j1_Start + 1);
        string idx_string = std::to_string(j1_Start) + "_" + std::to_string(j2_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int j1LB6 = j1_Start;
        for ( int j1 = j1LB6; j1 < 1023; j1++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrsymmat_addr25( j1) == calAddrsymmat_addr32(j1_Start, j2_Start)) {
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
                    if ( calAddrsymmat_addr26( j1, j2) == calAddrsymmat_addr32(j1_Start, j2_Start)) {
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
                        if ( calAddrsymmat_addr29( j1, j2, i) == calAddrsymmat_addr32(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsymmat_addr30( j1, j2, i) == calAddrsymmat_addr32(j1_Start, j2_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr31( j1, j2) == calAddrsymmat_addr32(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrsymmat_addr32( j1, j2) == calAddrsymmat_addr32(j1_Start, j2_Start)) {
                        rtHistoCal(cnt);
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
            if ( calAddrsymmat_addr33( ) == calAddrsymmat_addr32(j1_Start, j2_Start)) {
                rtHistoCal(cnt);
                goto EndSample;
            }
        }
EndSample:
        s++;
        }
}
void ref_symmat_addr33() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
    /* Generating reuse search code */

    s++;
    }
}
int main() {
    ref_stddev_addr6();
    ref_mean_addr0();
    ref_data_addr1();
    ref_mean_addr10();
    ref_stddev_addr11();
    ref_stddev_addr12();
    ref_mean_addr2();
    ref_mean_addr3();
    ref_mean_addr4();
    ref_mean_addr5();
    ref_stddev_addr17();
    ref_stddev_addr18();
    ref_data_addr7();
    ref_mean_addr8();
    ref_data_addr9();
    ref_symmat_addr26();
    ref_data_addr27();
    ref_data_addr28();
    ref_symmat_addr29();
    ref_symmat_addr30();
    ref_stddev_addr13();
    ref_stddev_addr14();
    ref_symmat_addr31();
    ref_stddev_addr15();
    ref_stddev_addr16();
    ref_mean_addr19();
    ref_data_addr20();
    ref_data_addr21();
    ref_stddev_addr22();
    ref_data_addr23();
    ref_data_addr24();
    ref_symmat_addr25();
    ref_symmat_addr32();
    ref_symmat_addr33();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: correlation */ 
