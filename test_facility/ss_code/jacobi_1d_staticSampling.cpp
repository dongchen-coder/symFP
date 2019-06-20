
 /* Start to analysis array index
Array index info
A.addr (i + 1)
A.addr (i - 1)
A.addr i
B.addr i
B.addr (i - 1)
B.addr i
B.addr (i + 1)
A.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %B

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (0, 10)
--Loop inc: (t + 1)
--Loop predicate: <
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------array access A.addr (i - 1)
------array access A.addr i
------array access A.addr (i + 1)
------array access B.addr i
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------array access B.addr (i - 1)
------array access B.addr i
------array access B.addr (i + 1)
------array access A.addr i

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 0
------Sample number: 1
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
int calAddrA_addr0( int t, int i) {
    int result = ((i - 1)) * 8 / 64;
    return result;
}
int calAddrA_addr1( int t, int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrA_addr2( int t, int i) {
    int result = ((i + 1)) * 8 / 64;
    return result;
}
int calAddrB_addr0( int t, int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrB_addr1( int t, int i) {
    int result = ((i - 1)) * 8 / 64;
    return result;
}
int calAddrB_addr2( int t, int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrB_addr3( int t, int i) {
    int result = ((i + 1)) * 8 / 64;
    return result;
}
int calAddrA_addr3( int t, int i) {
    int result = (i) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( t, i) == calAddrA_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( t, i) == calAddrA_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( t, i) == calAddrA_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( t, i) == calAddrA_addr0(t_Start, i_Start)) {
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
void ref_A_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( t, i) == calAddrA_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( t, i) == calAddrA_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( t, i) == calAddrA_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( t, i) == calAddrA_addr1(t_Start, i_Start)) {
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
void ref_A_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( t, i) == calAddrA_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( t, i) == calAddrA_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( t, i) == calAddrA_addr2(t_Start, i_Start)) {
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
            for ( int i = iLB2; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( t, i) == calAddrA_addr2(t_Start, i_Start)) {
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
void ref_A_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( t, i) == calAddrA_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( t, i) == calAddrA_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( t, i) == calAddrA_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB2 = 1;
            if ( t == t_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( t, i) == calAddrA_addr3(t_Start, i_Start)) {
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
void ref_B_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr0( t, i) == calAddrB_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            {
            int iLB2 = 1;
            for ( int i = iLB2; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr1( t, i) == calAddrB_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr2( t, i) == calAddrB_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr3( t, i) == calAddrB_addr0(t_Start, i_Start)) {
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
void ref_B_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr0( t, i) == calAddrB_addr1(t_Start, i_Start)) {
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
            for ( int i = iLB2; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr1( t, i) == calAddrB_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr2( t, i) == calAddrB_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr3( t, i) == calAddrB_addr1(t_Start, i_Start)) {
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
void ref_B_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr0( t, i) == calAddrB_addr2(t_Start, i_Start)) {
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
            for ( int i = iLB2; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr1( t, i) == calAddrB_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr2( t, i) == calAddrB_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr3( t, i) == calAddrB_addr2(t_Start, i_Start)) {
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
void ref_B_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t < 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr0( t, i) == calAddrB_addr3(t_Start, i_Start)) {
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
            for ( int i = iLB2; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr1( t, i) == calAddrB_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr2( t, i) == calAddrB_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrB_addr3( t, i) == calAddrB_addr3(t_Start, i_Start)) {
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
int main() {
    ref_A_addr0();
    ref_A_addr1();
    ref_A_addr2();
    ref_A_addr3();
    ref_B_addr0();
    ref_B_addr1();
    ref_B_addr2();
    ref_B_addr3();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: jacobi_1d */ 
