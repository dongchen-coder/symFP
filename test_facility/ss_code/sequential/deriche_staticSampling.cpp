
 /* Start to analysis array index
Array index info: Total number of references: 20
y1.addr ((i * 1024) + j)
imgIn.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
imgIn.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgIn.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %y1
double* %imgIn
double* %y2
double* %imgOut
double %alpha

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
------array access imgIn.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
------array access imgIn.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (1023, 0)
----Loop inc: (j + -1)
----Loop predicate: >=
------array access y2.addr ((i * 1024) + j)
------array access imgIn.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access y1.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)
--j
--Loop Bound: (0, 1024)
--Loop inc: (j + 1)
--Loop predicate: <
----i
----Loop Bound: (0, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------array access imgOut.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
--j
--Loop Bound: (0, 1024)
--Loop inc: (j + 1)
--Loop predicate: <
----i
----Loop Bound: (1023, 0)
----Loop inc: (i + -1)
----Loop predicate: >=
------array access y2.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access y1.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
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
/* imgIn_addr ((i * 1024) + j) 0 */
int calAddrimgIn_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 1 */
int calAddry1_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgIn_addr ((i * 1024) + j) 2 */
int calAddrimgIn_addr2( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 3 */
int calAddry1_addr3( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 4 */
int calAddry2_addr4( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgIn_addr ((i * 1024) + j) 5 */
int calAddrimgIn_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 6 */
int calAddry2_addr6( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 7 */
int calAddry1_addr7( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 8 */
int calAddry2_addr8( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 9 */
int calAddrimgOut_addr9( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 10 */
int calAddrimgOut_addr10( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 11 */
int calAddry1_addr11( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 12 */
int calAddrimgOut_addr12( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 13 */
int calAddry1_addr13( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 14 */
int calAddry2_addr14( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 15 */
int calAddrimgOut_addr15( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 16 */
int calAddry2_addr16( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 17 */
int calAddry1_addr17( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 18 */
int calAddry2_addr18( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 19 */
int calAddrimgOut_addr19( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_y1_addr1() {
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
                    if ( calAddry1_addr1( i, j) == calAddry1_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr3( i, j) == calAddry1_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
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
            int jLB3 = 1023;
            for ( int j = jLB3; j >= 0; j--) {
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
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr7( i, j) == calAddry1_addr1(i_Start, j_Start)) {
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr1(i_Start, j_Start)) {
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
void ref_imgIn_addr2() {
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
                    if ( calAddrimgIn_addr0( i, j) == calAddrimgIn_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr2( i, j) == calAddrimgIn_addr2(i_Start, j_Start)) {
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
        int iLB2 = 0;
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 1023;
            for ( int j = jLB3; j >= 0; j--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr5( i, j) == calAddrimgIn_addr2(i_Start, j_Start)) {
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
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
void ref_y1_addr3() {
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
                    if ( calAddry1_addr1( i, j) == calAddry1_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr3( i, j) == calAddry1_addr3(i_Start, j_Start)) {
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
        int iLB2 = 0;
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 1023;
            for ( int j = jLB3; j >= 0; j--) {
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
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr7( i, j) == calAddry1_addr3(i_Start, j_Start)) {
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr3(i_Start, j_Start)) {
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
void ref_y2_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0 + 1) + 0;
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
            int jLB3 = 1023;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j >= 0; j--) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr4( i, j) == calAddry2_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr6( i, j) == calAddry2_addr6(i_Start, j_Start)) {
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
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr8( i, j) == calAddry2_addr6(i_Start, j_Start)) {
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr6(i_Start, j_Start)) {
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
void ref_y1_addr7() {
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
                    if ( calAddry1_addr7( i, j) == calAddry1_addr7(i_Start, j_Start)) {
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr7(i_Start, j_Start)) {
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
void ref_imgIn_addr0() {
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
                    if ( calAddrimgIn_addr0( i, j) == calAddrimgIn_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr2( i, j) == calAddrimgIn_addr0(i_Start, j_Start)) {
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
        int iLB2 = 0;
        for ( int i = iLB2; i < 1024; i++) {
            {
            int jLB3 = 1023;
            for ( int j = jLB3; j >= 0; j--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr5( i, j) == calAddrimgIn_addr0(i_Start, j_Start)) {
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
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
void ref_y2_addr8() {
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
                    if ( calAddry2_addr8( i, j) == calAddry2_addr8(i_Start, j_Start)) {
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr8(i_Start, j_Start)) {
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
void ref_imgOut_addr9() {
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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr9( i, j) == calAddrimgOut_addr9(i_Start, j_Start)) {
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr10( j, i) == calAddrimgOut_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr12( j, i) == calAddrimgOut_addr9(i_Start, j_Start)) {
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
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr15( j, i) == calAddrimgOut_addr9(i_Start, j_Start)) {
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
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr9(i_Start, j_Start)) {
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
void ref_y1_addr11() {
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
        int jLB6 = j_Start;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr11(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr11(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr11(j_Start, i_Start)) {
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
void ref_imgOut_addr12() {
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
        int jLB6 = j_Start;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr10( j, i) == calAddrimgOut_addr12(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr12( j, i) == calAddrimgOut_addr12(j_Start, i_Start)) {
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
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr15( j, i) == calAddrimgOut_addr12(j_Start, i_Start)) {
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
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr12(j_Start, i_Start)) {
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
void ref_y1_addr13() {
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
        int jLB6 = j_Start;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr11( j, i) == calAddry1_addr13(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr13( j, i) == calAddry1_addr13(j_Start, i_Start)) {
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
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr13(j_Start, i_Start)) {
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
void ref_y2_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0 + 1) + 0;
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
            int jLB3 = 1023;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j >= 0; j--) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr4( i, j) == calAddry2_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr6( i, j) == calAddry2_addr4(i_Start, j_Start)) {
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
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr8( i, j) == calAddry2_addr4(i_Start, j_Start)) {
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
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr4(i_Start, j_Start)) {
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
void ref_imgIn_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0 + 1) + 0;
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
            int jLB3 = 1023;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j >= 0; j--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgIn_addr5( i, j) == calAddrimgIn_addr5(i_Start, j_Start)) {
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
            int jLB5 = 0;
            for ( int j = jLB5; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB6 = 0;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
void ref_imgOut_addr10() {
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
        int jLB6 = j_Start;
        for ( int j = jLB6; j < 1024; j++) {
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 1024; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr10( j, i) == calAddrimgOut_addr10(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr12( j, i) == calAddrimgOut_addr10(j_Start, i_Start)) {
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
        int jLB8 = 0;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr15( j, i) == calAddrimgOut_addr10(j_Start, i_Start)) {
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
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr10(j_Start, i_Start)) {
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
void ref_y2_addr14() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB8 = j_Start;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            if ( j == j_Start ) {
                iLB9 = i_Start;
            }
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr14(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr14(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
        }
        }
        {
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr14(j_Start, i_Start)) {
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
void ref_imgOut_addr15() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB8 = j_Start;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            if ( j == j_Start ) {
                iLB9 = i_Start;
            }
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr15( j, i) == calAddrimgOut_addr15(j_Start, i_Start)) {
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
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr15(j_Start, i_Start)) {
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
void ref_y2_addr16() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int jLB8 = j_Start;
        for ( int j = jLB8; j < 1024; j++) {
            {
            int iLB9 = 1023;
            if ( j == j_Start ) {
                iLB9 = i_Start;
            }
            for ( int i = iLB9; i >= 0; i--) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr14( j, i) == calAddry2_addr16(j_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr16( j, i) == calAddry2_addr16(j_Start, i_Start)) {
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
        int iLB10 = 0;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr16(j_Start, i_Start)) {
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
void ref_y1_addr17() {
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
        int iLB10 = i_Start;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            if ( i == i_Start ) {
                jLB11 = j_Start;
            }
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr17( i, j) == calAddry1_addr17(i_Start, j_Start)) {
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
EndSample:
        s++;
        }
}
void ref_y2_addr18() {
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
        int iLB10 = i_Start;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            if ( i == i_Start ) {
                jLB11 = j_Start;
            }
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr18( i, j) == calAddry2_addr18(i_Start, j_Start)) {
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
EndSample:
        s++;
        }
}
void ref_imgOut_addr19() {
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
        int iLB10 = i_Start;
        for ( int i = iLB10; i < 1024; i++) {
            {
            int jLB11 = 0;
            if ( i == i_Start ) {
                jLB11 = j_Start;
            }
            for ( int j = jLB11; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr19( i, j) == calAddrimgOut_addr19(i_Start, j_Start)) {
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
    ref_y1_addr1();
    ref_imgIn_addr2();
    ref_y1_addr3();
    ref_y2_addr6();
    ref_y1_addr7();
    ref_imgIn_addr0();
    ref_y2_addr8();
    ref_imgOut_addr9();
    ref_y1_addr11();
    ref_imgOut_addr12();
    ref_y1_addr13();
    ref_y2_addr4();
    ref_imgIn_addr5();
    ref_imgOut_addr10();
    ref_y2_addr14();
    ref_imgOut_addr15();
    ref_y2_addr16();
    ref_y1_addr17();
    ref_y2_addr18();
    ref_imgOut_addr19();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: deriche */ 
