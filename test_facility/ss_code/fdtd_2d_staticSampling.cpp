
 /* Start to analysis array index
Array index info
ey.addr (0 + j)
_fict_.addr t
ey.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)
hz.addr (((i - 1) * 1024) + j)
ey.addr ((i * 1024) + j)
ex.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)
hz.addr (((i * 1024) + j) - 1)
ex.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)
ex.addr (((i * 1024) + j) + 1)
ex.addr ((i * 1024) + j)
ey.addr (((i + 1) * 1024) + j)
ey.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %_fict_
double* %ey
double* %ex
double* %hz

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (0, 10)
--Loop inc: (t + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access _fict_.addr t
------array access ey.addr (0 + j)
----i
----Loop Bound: (1, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (0, 1024)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access ey.addr ((i * 1024) + j)
--------array access hz.addr ((i * 1024) + j)
--------array access hz.addr (((i - 1) * 1024) + j)
--------array access ey.addr ((i * 1024) + j)
----i
----Loop Bound: (0, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (1, 1024)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access ex.addr ((i * 1024) + j)
--------array access hz.addr ((i * 1024) + j)
--------array access hz.addr (((i * 1024) + j) - 1)
--------array access ex.addr ((i * 1024) + j)
----i
----Loop Bound: (0, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (0, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access hz.addr ((i * 1024) + j)
--------array access ex.addr (((i * 1024) + j) + 1)
--------array access ex.addr ((i * 1024) + j)
--------array access ey.addr (((i + 1) * 1024) + j)
--------array access ey.addr ((i * 1024) + j)
--------array access hz.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 0
------Sample number: 1
------Sample number: 1
--------Sample number: 10
------Sample number: 1
--------Sample number: 10
------Sample number: 1
--------Sample number: 10
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
int calAddr_fict__addr0( int t, int j) {
    int result = (t) * 8 / 64;
    return result;
}
int calAddrey_addr0( int t, int j) {
    int result = ((0 + j)) * 8 / 64;
    return result;
}
int calAddrey_addr1( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrhz_addr0( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrhz_addr1( int t, int i, int j) {
    int result = ((((i - 1) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrey_addr2( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrex_addr0( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrhz_addr2( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrhz_addr3( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrex_addr1( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrhz_addr4( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrex_addr2( int t, int i, int j) {
    int result = ((((i * 1024) + j) + 1)) * 8 / 64;
    return result;
}
int calAddrex_addr3( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrey_addr3( int t, int i, int j) {
    int result = ((((i + 1) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrey_addr4( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrhz_addr5( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref__fict__addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            if ( t == t_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddr_fict__addr0( t, j) == calAddr_fict__addr0(t_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
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
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
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
void ref_ex_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
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
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1024; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr0( t, i, j) == calAddrex_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr1( t, i, j) == calAddrex_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr2( t, i, j) == calAddrex_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr3( t, i, j) == calAddrex_addr0(t_Start, i_Start, j_Start)) {
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
void ref_ex_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
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
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1024; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr0( t, i, j) == calAddrex_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr1( t, i, j) == calAddrex_addr1(t_Start, i_Start, j_Start)) {
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
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr2( t, i, j) == calAddrex_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr3( t, i, j) == calAddrex_addr1(t_Start, i_Start, j_Start)) {
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
void ref_ex_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
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
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr0( t, i, j) == calAddrex_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr1( t, i, j) == calAddrex_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr2( t, i, j) == calAddrex_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr3( t, i, j) == calAddrex_addr2(t_Start, i_Start, j_Start)) {
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
void ref_ex_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
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
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr0( t, i, j) == calAddrex_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr1( t, i, j) == calAddrex_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr2( t, i, j) == calAddrex_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr3( t, i, j) == calAddrex_addr3(t_Start, i_Start, j_Start)) {
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
void ref_ey_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            if ( t == t_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( t, j) == calAddrey_addr0(t_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( t, i, j) == calAddrey_addr0(t_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( t, i, j) == calAddrey_addr0(t_Start, j_Start)) {
                            rtHistoCal(cnt);
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
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( t, i, j) == calAddrey_addr0(t_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( t, i, j) == calAddrey_addr0(t_Start, j_Start)) {
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
void ref_ey_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 1) + 1;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( t, j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int iLB2 = 1;
            if ( t == t_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( t, i, j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( t, i, j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( t, i, j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( t, i, j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
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
void ref_ey_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 1) + 1;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( t, j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int iLB2 = 1;
            if ( t == t_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( t, i, j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( t, i, j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
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
            int iLB4 = 0;
            for ( int i = iLB4; i < 1024; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( t, i, j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( t, i, j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
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
void ref_ey_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( t, j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( t, i, j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( t, i, j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( t, i, j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( t, i, j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
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
void ref_ey_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( t, j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( t, i, j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( t, i, j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( t, i, j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( t, i, j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
void ref_hz_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 1) + 1;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
            if ( t == t_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( t, i, j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( t, i, j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 0;
            for ( int i = iLB4; i < 1024; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( t, i, j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( t, i, j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( t, i, j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( t, i, j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
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
void ref_hz_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 1) + 1;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
            if ( t == t_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( t, i, j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( t, i, j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( t, i, j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( t, i, j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( t, i, j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( t, i, j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
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
void ref_hz_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( t, i, j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( t, i, j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 0;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1024; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( t, i, j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( t, i, j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( t, i, j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( t, i, j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
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
void ref_hz_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( t, i, j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( t, i, j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 0;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1024; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( t, i, j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( t, i, j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
            int iLB6 = 0;
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( t, i, j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( t, i, j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
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
void ref_hz_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( t, i, j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( t, i, j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 0;
            for ( int i = iLB4; i < 1024; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( t, i, j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( t, i, j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( t, i, j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
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
                        if ( calAddrhz_addr5( t, i, j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
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
void ref_hz_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1024; i++) {
                {
                int jLB3 = 0;
                for ( int j = jLB3; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( t, i, j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( t, i, j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 0;
            for ( int i = iLB4; i < 1024; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1024; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( t, i, j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( t, i, j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( t, i, j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( t, i, j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
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
    ref__fict__addr0();
    ref_ex_addr0();
    ref_ex_addr1();
    ref_ex_addr2();
    ref_ex_addr3();
    ref_ey_addr0();
    ref_ey_addr1();
    ref_ey_addr2();
    ref_ey_addr3();
    ref_ey_addr4();
    ref_hz_addr0();
    ref_hz_addr1();
    ref_hz_addr2();
    ref_hz_addr3();
    ref_hz_addr4();
    ref_hz_addr5();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: fdtd_2d */ 
