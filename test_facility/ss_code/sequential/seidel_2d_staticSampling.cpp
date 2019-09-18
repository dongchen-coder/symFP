
 /* Start to analysis array index
Array index info: Total number of references: 10
A.addr ((((i - 1) * 1024) + j) - 1)
A.addr ((i * 1024) + j)
A.addr (((i * 1024) + j) + 1)
A.addr (((i - 1) * 1024) + j)
A.addr (((i + 1) * 1024) + j)
A.addr ((((i + 1) * 1024) + j) + 1)
A.addr ((((i - 1) * 1024) + j) + 1)
A.addr (((i * 1024) + j) - 1)
A.addr ((((i + 1) * 1024) + j) - 1)
A.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (0, 9)
--Loop inc: (t + 1)
--Loop predicate: <=
----i
----Loop Bound: (1, 1022)
----Loop inc: (i + 1)
----Loop predicate: <=
------j
------Loop Bound: (1, 1022)
------Loop inc: (j + 1)
------Loop predicate: <=
--------array access A.addr ((((i - 1) * 1024) + j) - 1)
--------array access A.addr (((i - 1) * 1024) + j)
--------array access A.addr ((((i - 1) * 1024) + j) + 1)
--------array access A.addr (((i * 1024) + j) - 1)
--------array access A.addr ((i * 1024) + j)
--------array access A.addr (((i * 1024) + j) + 1)
--------array access A.addr ((((i + 1) * 1024) + j) - 1)
--------array access A.addr (((i + 1) * 1024) + j)
--------array access A.addr ((((i + 1) * 1024) + j) + 1)
--------array access A.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 0
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
/* A_addr ((((i - 1) * 1024) + j) - 1) 0 */
int calAddrA_addr0( int t, int i, int j) {
    int result = (((((i - 1) * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* A_addr (((i - 1) * 1024) + j) 1 */
int calAddrA_addr1( int t, int i, int j) {
    int result = ((((i - 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr ((((i - 1) * 1024) + j) + 1) 2 */
int calAddrA_addr2( int t, int i, int j) {
    int result = (((((i - 1) * 1024) + j) + 1)) * 8 / 64;
    return result;
}
/* A_addr (((i * 1024) + j) - 1) 3 */
int calAddrA_addr3( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + j) 4 */
int calAddrA_addr4( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr (((i * 1024) + j) + 1) 5 */
int calAddrA_addr5( int t, int i, int j) {
    int result = ((((i * 1024) + j) + 1)) * 8 / 64;
    return result;
}
/* A_addr ((((i + 1) * 1024) + j) - 1) 6 */
int calAddrA_addr6( int t, int i, int j) {
    int result = (((((i + 1) * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* A_addr (((i + 1) * 1024) + j) 7 */
int calAddrA_addr7( int t, int i, int j) {
    int result = ((((i + 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr ((((i + 1) * 1024) + j) + 1) 8 */
int calAddrA_addr8( int t, int i, int j) {
    int result = (((((i + 1) * 1024) + j) + 1)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + j) 9 */
int calAddrA_addr9( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr0(t_Start, i_Start, j_Start)) {
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
void ref_A_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr4(t_Start, i_Start, j_Start)) {
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
void ref_A_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr5(t_Start, i_Start, j_Start)) {
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
void ref_A_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr1(t_Start, i_Start, j_Start)) {
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
void ref_A_addr7() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr7(t_Start, i_Start, j_Start)) {
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
void ref_A_addr8() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr8(t_Start, i_Start, j_Start)) {
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
void ref_A_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr2(t_Start, i_Start, j_Start)) {
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
void ref_A_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr3(t_Start, i_Start, j_Start)) {
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
void ref_A_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr6(t_Start, i_Start, j_Start)) {
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
void ref_A_addr9() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (9 - 0 + 1) + 0;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1022 - 1 + 1) + 1;
        if ( (1022 - 1 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 9; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i <= 1022; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j <= 1022; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr4( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr5( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr6( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr7( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr8( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr9( t, i, j) == calAddrA_addr9(t_Start, i_Start, j_Start)) {
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
    ref_A_addr0();
    ref_A_addr4();
    ref_A_addr5();
    ref_A_addr1();
    ref_A_addr7();
    ref_A_addr8();
    ref_A_addr2();
    ref_A_addr3();
    ref_A_addr6();
    ref_A_addr9();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: seidel_2d */ 
