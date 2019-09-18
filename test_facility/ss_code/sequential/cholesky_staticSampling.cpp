
 /* Start to analysis array index
Array index info: Total number of references: 13
A.addr ((i * 1024) + k)
A.addr ((j * 1024) + k)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + i)
A.addr ((j * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + i)
A.addr ((i * 1024) + i)
A.addr ((i * 1024) + i)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, i)
----Loop inc: (j + 1)
----Loop predicate: <
------k
------Loop Bound: (0, j)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 1024) + k)
--------array access A.addr ((j * 1024) + k)
--------array access A.addr ((i * 1024) + j)
--------array access A.addr ((i * 1024) + j)
------array access A.addr ((j * 1024) + j)
------array access A.addr ((i * 1024) + j)
------array access A.addr ((i * 1024) + j)
----k
----Loop Bound: (0, i)
----Loop inc: (k + 1)
----Loop predicate: <
------array access A.addr ((i * 1024) + k)
------array access A.addr ((i * 1024) + k)
------array access A.addr ((i * 1024) + i)
------array access A.addr ((i * 1024) + i)
----array access A.addr ((i * 1024) + i)
----array access A.addr ((i * 1024) + i)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 0 0 
Dump stride: 1 1 
init counter: 0 0 0 
Dump stride: 1 1 1 
init counter: 0 0 
Dump stride: 1 1 
Dump tree:
----Sample number: 10
------Sample number: 52
--------Sample number: 178
------Sample number: 52
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
/* A_addr ((i * 1024) + k) 0 */
int calAddrA_addr0( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* A_addr ((j * 1024) + k) 1 */
int calAddrA_addr1( int i, int j, int k) {
    int result = (((j * 1024) + k)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + j) 2 */
int calAddrA_addr2( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + j) 3 */
int calAddrA_addr3( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr ((j * 1024) + j) 4 */
int calAddrA_addr4( int i, int j) {
    int result = (((j * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + j) 5 */
int calAddrA_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + j) 6 */
int calAddrA_addr6( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + k) 7 */
int calAddrA_addr7( int i, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + k) 8 */
int calAddrA_addr8( int i, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + i) 9 */
int calAddrA_addr9( int i, int k) {
    int result = (((i * 1024) + i)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + i) 10 */
int calAddrA_addr10( int i, int k) {
    int result = (((i * 1024) + i)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + i) 11 */
int calAddrA_addr11( int i) {
    int result = (((i * 1024) + i)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + i) 12 */
int calAddrA_addr12( int i) {
    int result = (((i * 1024) + i)) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        if ( (j_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
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
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        if ( (j_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
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
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        if ( (j_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
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
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        if ( (j_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr11() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr11(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr11(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr11(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr11(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr11(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr11(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr4(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr4(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr5(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr5(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr6(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr6(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr7() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr7(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr7(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr7(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr7(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            if ( i == i_Start ) {
                kLB3 = k_Start;
            }
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr7(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr7(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr8() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr8(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr8(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr8(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr8(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            if ( i == i_Start ) {
                kLB3 = k_Start;
            }
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr8(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr8(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr9() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr9(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr9(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr9(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr9(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            if ( i == i_Start ) {
                kLB3 = k_Start;
            }
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr9(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr9(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr10() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr10(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr10(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr10(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr10(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            if ( i == i_Start ) {
                kLB3 = k_Start;
            }
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr10(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr10(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr12() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr12(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr12(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr12(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr12(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr12(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr12(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
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
    ref_A_addr11();
    ref_A_addr4();
    ref_A_addr5();
    ref_A_addr6();
    ref_A_addr7();
    ref_A_addr8();
    ref_A_addr9();
    ref_A_addr10();
    ref_A_addr12();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: cholesky */ 
