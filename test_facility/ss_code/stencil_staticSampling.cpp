
 /* Start to analysis array index
Array index info
a.addr (((1026 * i) + j) + 1)
a.addr ((1026 * i) + j)
a.addr (((1026 * i) + j) - 1)
a.addr ((1026 * (i - 1)) + j)
a.addr ((1026 * (i + 1)) + j)
b.addr ((1026 * i) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %a
double* %b

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (1, 1025)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (1, 1025)
----Loop inc: (j + 1)
----Loop predicate: <
------array access a.addr ((1026 * i) + j)
------array access a.addr (((1026 * i) + j) + 1)
------array access a.addr (((1026 * i) + j) - 1)
------array access a.addr ((1026 * (i - 1)) + j)
------array access a.addr ((1026 * (i + 1)) + j)
------array access b.addr ((1026 * i) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
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
int calAddra_addr0( int i, int j) {
    int result = (((1026 * i) + j)) * 8 / 64;
    return result;
}
int calAddra_addr1( int i, int j) {
    int result = ((((1026 * i) + j) + 1)) * 8 / 64;
    return result;
}
int calAddra_addr2( int i, int j) {
    int result = ((((1026 * i) + j) - 1)) * 8 / 64;
    return result;
}
int calAddra_addr3( int i, int j) {
    int result = (((1026 * (i - 1)) + j)) * 8 / 64;
    return result;
}
int calAddra_addr4( int i, int j) {
    int result = (((1026 * (i + 1)) + j)) * 8 / 64;
    return result;
}
int calAddrb_addr0( int i, int j) {
    int result = (((1026 * i) + j)) * 8 / 64;
    return result;
}
void ref_a_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1025; i++) {
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1025; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr0( i, j) == calAddra_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr0(i_Start, j_Start)) {
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
void ref_a_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1025; i++) {
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1025; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr0( i, j) == calAddra_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr1(i_Start, j_Start)) {
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
void ref_a_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1025; i++) {
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1025; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr0( i, j) == calAddra_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr2(i_Start, j_Start)) {
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
void ref_a_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1025; i++) {
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1025; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr0( i, j) == calAddra_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr3(i_Start, j_Start)) {
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
void ref_a_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1025; i++) {
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1025; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr0( i, j) == calAddra_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr4(i_Start, j_Start)) {
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
void ref_b_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1025; i++) {
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1025; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrb_addr0( i, j) == calAddrb_addr0(i_Start, j_Start)) {
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
    ref_a_addr0();
    ref_a_addr1();
    ref_a_addr2();
    ref_a_addr3();
    ref_a_addr4();
    ref_b_addr0();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: stencil */ 
