
 /* Start to analysis array index
Array index info
b.addr i
L.addr ((i * 1024) + j)
x.addr j
x.addr i
x.addr i
x.addr i
x.addr i
L.addr ((i * 1024) + i)
x.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %x
double* %b
double* %L

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----array access b.addr i
----array access x.addr i
----j
----Loop Bound: (0, i)
----Loop inc: (j + 1)
----Loop predicate: <
------array access L.addr ((i * 1024) + j)
------array access x.addr j
------array access x.addr i
------array access x.addr i
----array access x.addr i
----array access L.addr ((i * 1024) + i)
----array access x.addr i

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 0 0 
Dump stride: 1 1 
Dump tree:
----Sample number: 10
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
int calAddrb_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrx_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrL_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrx_addr1( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrx_addr2( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrx_addr3( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrx_addr4( int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrL_addr1( int i) {
    int result = (((i * 1024) + i)) * 8 / 64;
    return result;
}
int calAddrx_addr5( int i) {
    int result = (i) * 8 / 64;
    return result;
}
void ref_L_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_L_addr0 = -1;
    uint64_t prev_i_Start_L_addr0 = -1;
    uint64_t prev_i_End_L_addr0 = -1;
    uint64_t prev_j_Start_L_addr0 = -1;
    uint64_t prev_j_End_L_addr0 = -1;
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

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_L_addr0 != -1) {
            if ( calAddrL_addr0( i_Start - prev_i_Start_L_addr0 + prev_i_End_L_addr0, j_Start - prev_j_Start_L_addr0 + prev_j_End_L_addr0) == calAddrL_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_L_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrL_addr0( i, j) == calAddrL_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_L_addr0 = cnt;
                        prev_i_Start_L_addr0 = i_Start;
                        prev_i_End_L_addr0 = i;
                        prev_j_Start_L_addr0 = j_Start;
                        prev_j_End_L_addr0 = j;
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
            if (cntStart == true) {
                cnt++;
                if ( calAddrL_addr1( i) == calAddrL_addr0(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
        }
        }
EndSample:
        s++;
        }
}
void ref_L_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_L_addr1 = -1;
    uint64_t prev_i_Start_L_addr1 = -1;
    uint64_t prev_i_End_L_addr1 = -1;
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

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_L_addr1 != -1) {
            if ( calAddrL_addr1( i_Start - prev_i_Start_L_addr1 + prev_i_End_L_addr1) == calAddrL_addr1(i_Start)) {
                rtHistoCal(prev_cnt_L_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrL_addr0( i, j) == calAddrL_addr1(i_Start)) {
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
            if (cntStart == true) {
                cnt++;
                if ( calAddrL_addr1( i) == calAddrL_addr1(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_L_addr1 = cnt;
                    prev_i_Start_L_addr1 = i_Start;
                    prev_i_End_L_addr1 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
        }
        }
EndSample:
        s++;
        }
}
void ref_b_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_b_addr0 = -1;
    uint64_t prev_i_Start_b_addr0 = -1;
    uint64_t prev_i_End_b_addr0 = -1;
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

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_b_addr0 != -1) {
            if ( calAddrb_addr0( i_Start - prev_i_Start_b_addr0 + prev_i_End_b_addr0) == calAddrb_addr0(i_Start)) {
                rtHistoCal(prev_cnt_b_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrb_addr0( i) == calAddrb_addr0(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_b_addr0 = cnt;
                    prev_i_Start_b_addr0 = i_Start;
                    prev_i_End_b_addr0 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
        }
        }
EndSample:
        s++;
        }
}
void ref_x_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr0 = -1;
    uint64_t prev_i_Start_x_addr0 = -1;
    uint64_t prev_i_End_x_addr0 = -1;
    uint64_t prev_cnt_x_addr4 = -1;
    uint64_t prev_i_Start_x_addr4 = -1;
    uint64_t prev_i_End_x_addr4 = -1;
    uint64_t prev_cnt_x_addr5 = -1;
    uint64_t prev_i_Start_x_addr5 = -1;
    uint64_t prev_i_End_x_addr5 = -1;
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

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_x_addr0 != -1) {
            if ( calAddrx_addr0( i_Start - prev_i_Start_x_addr0 + prev_i_End_x_addr0) == calAddrx_addr0(i_Start)) {
                rtHistoCal(prev_cnt_x_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr4 != -1) {
            if ( calAddrx_addr4( i_Start - prev_i_Start_x_addr4 + prev_i_End_x_addr4) == calAddrx_addr0(i_Start)) {
                rtHistoCal(prev_cnt_x_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr5 != -1) {
            if ( calAddrx_addr5( i_Start - prev_i_Start_x_addr5 + prev_i_End_x_addr5) == calAddrx_addr0(i_Start)) {
                rtHistoCal(prev_cnt_x_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr0( i) == calAddrx_addr0(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr0 = cnt;
                    prev_i_Start_x_addr0 = i_Start;
                    prev_i_End_x_addr0 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr0(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr2( i, j) == calAddrx_addr0(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr3( i, j) == calAddrx_addr0(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr4( i) == calAddrx_addr0(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr4 = cnt;
                    prev_i_Start_x_addr4 = i_Start;
                    prev_i_End_x_addr4 = i;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr5( i) == calAddrx_addr0(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr5 = cnt;
                    prev_i_Start_x_addr5 = i_Start;
                    prev_i_End_x_addr5 = i;
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_x_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr1 = -1;
    uint64_t prev_i_Start_x_addr1 = -1;
    uint64_t prev_i_End_x_addr1 = -1;
    uint64_t prev_j_Start_x_addr1 = -1;
    uint64_t prev_j_End_x_addr1 = -1;
    uint64_t prev_cnt_x_addr2 = -1;
    uint64_t prev_i_Start_x_addr2 = -1;
    uint64_t prev_i_End_x_addr2 = -1;
    uint64_t prev_j_Start_x_addr2 = -1;
    uint64_t prev_j_End_x_addr2 = -1;
    uint64_t prev_cnt_x_addr3 = -1;
    uint64_t prev_i_Start_x_addr3 = -1;
    uint64_t prev_i_End_x_addr3 = -1;
    uint64_t prev_j_Start_x_addr3 = -1;
    uint64_t prev_j_End_x_addr3 = -1;
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

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_x_addr1 != -1) {
            if ( calAddrx_addr1( i_Start - prev_i_Start_x_addr1 + prev_i_End_x_addr1, j_Start - prev_j_Start_x_addr1 + prev_j_End_x_addr1) == calAddrx_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr2 != -1) {
            if ( calAddrx_addr2( i_Start - prev_i_Start_x_addr2 + prev_i_End_x_addr2, j_Start - prev_j_Start_x_addr2 + prev_j_End_x_addr2) == calAddrx_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr3 != -1) {
            if ( calAddrx_addr3( i_Start - prev_i_Start_x_addr3 + prev_i_End_x_addr3, j_Start - prev_j_Start_x_addr3 + prev_j_End_x_addr3) == calAddrx_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr0( i) == calAddrx_addr1(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr1 = cnt;
                        prev_i_Start_x_addr1 = i_Start;
                        prev_i_End_x_addr1 = i;
                        prev_j_Start_x_addr1 = j_Start;
                        prev_j_End_x_addr1 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr2( i, j) == calAddrx_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr2 = cnt;
                        prev_i_Start_x_addr2 = i_Start;
                        prev_i_End_x_addr2 = i;
                        prev_j_Start_x_addr2 = j_Start;
                        prev_j_End_x_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr3( i, j) == calAddrx_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr3 = cnt;
                        prev_i_Start_x_addr3 = i_Start;
                        prev_i_End_x_addr3 = i;
                        prev_j_Start_x_addr3 = j_Start;
                        prev_j_End_x_addr3 = j;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr4( i) == calAddrx_addr1(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr5( i) == calAddrx_addr1(i_Start, j_Start)) {
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
void ref_x_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr1 = -1;
    uint64_t prev_i_Start_x_addr1 = -1;
    uint64_t prev_i_End_x_addr1 = -1;
    uint64_t prev_j_Start_x_addr1 = -1;
    uint64_t prev_j_End_x_addr1 = -1;
    uint64_t prev_cnt_x_addr2 = -1;
    uint64_t prev_i_Start_x_addr2 = -1;
    uint64_t prev_i_End_x_addr2 = -1;
    uint64_t prev_j_Start_x_addr2 = -1;
    uint64_t prev_j_End_x_addr2 = -1;
    uint64_t prev_cnt_x_addr3 = -1;
    uint64_t prev_i_Start_x_addr3 = -1;
    uint64_t prev_i_End_x_addr3 = -1;
    uint64_t prev_j_Start_x_addr3 = -1;
    uint64_t prev_j_End_x_addr3 = -1;
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

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_x_addr1 != -1) {
            if ( calAddrx_addr1( i_Start - prev_i_Start_x_addr1 + prev_i_End_x_addr1, j_Start - prev_j_Start_x_addr1 + prev_j_End_x_addr1) == calAddrx_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr2 != -1) {
            if ( calAddrx_addr2( i_Start - prev_i_Start_x_addr2 + prev_i_End_x_addr2, j_Start - prev_j_Start_x_addr2 + prev_j_End_x_addr2) == calAddrx_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr3 != -1) {
            if ( calAddrx_addr3( i_Start - prev_i_Start_x_addr3 + prev_i_End_x_addr3, j_Start - prev_j_Start_x_addr3 + prev_j_End_x_addr3) == calAddrx_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr0( i) == calAddrx_addr2(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr1 = cnt;
                        prev_i_Start_x_addr1 = i_Start;
                        prev_i_End_x_addr1 = i;
                        prev_j_Start_x_addr1 = j_Start;
                        prev_j_End_x_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr2( i, j) == calAddrx_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr2 = cnt;
                        prev_i_Start_x_addr2 = i_Start;
                        prev_i_End_x_addr2 = i;
                        prev_j_Start_x_addr2 = j_Start;
                        prev_j_End_x_addr2 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr3( i, j) == calAddrx_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr3 = cnt;
                        prev_i_Start_x_addr3 = i_Start;
                        prev_i_End_x_addr3 = i;
                        prev_j_Start_x_addr3 = j_Start;
                        prev_j_End_x_addr3 = j;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr4( i) == calAddrx_addr2(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr5( i) == calAddrx_addr2(i_Start, j_Start)) {
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
void ref_x_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr1 = -1;
    uint64_t prev_i_Start_x_addr1 = -1;
    uint64_t prev_i_End_x_addr1 = -1;
    uint64_t prev_j_Start_x_addr1 = -1;
    uint64_t prev_j_End_x_addr1 = -1;
    uint64_t prev_cnt_x_addr2 = -1;
    uint64_t prev_i_Start_x_addr2 = -1;
    uint64_t prev_i_End_x_addr2 = -1;
    uint64_t prev_j_Start_x_addr2 = -1;
    uint64_t prev_j_End_x_addr2 = -1;
    uint64_t prev_cnt_x_addr3 = -1;
    uint64_t prev_i_Start_x_addr3 = -1;
    uint64_t prev_i_End_x_addr3 = -1;
    uint64_t prev_j_Start_x_addr3 = -1;
    uint64_t prev_j_End_x_addr3 = -1;
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

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_x_addr1 != -1) {
            if ( calAddrx_addr1( i_Start - prev_i_Start_x_addr1 + prev_i_End_x_addr1, j_Start - prev_j_Start_x_addr1 + prev_j_End_x_addr1) == calAddrx_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr2 != -1) {
            if ( calAddrx_addr2( i_Start - prev_i_Start_x_addr2 + prev_i_End_x_addr2, j_Start - prev_j_Start_x_addr2 + prev_j_End_x_addr2) == calAddrx_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr3 != -1) {
            if ( calAddrx_addr3( i_Start - prev_i_Start_x_addr3 + prev_i_End_x_addr3, j_Start - prev_j_Start_x_addr3 + prev_j_End_x_addr3) == calAddrx_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_x_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr0( i) == calAddrx_addr3(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr1 = cnt;
                        prev_i_Start_x_addr1 = i_Start;
                        prev_i_End_x_addr1 = i;
                        prev_j_Start_x_addr1 = j_Start;
                        prev_j_End_x_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr2( i, j) == calAddrx_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr2 = cnt;
                        prev_i_Start_x_addr2 = i_Start;
                        prev_i_End_x_addr2 = i;
                        prev_j_Start_x_addr2 = j_Start;
                        prev_j_End_x_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr3( i, j) == calAddrx_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_x_addr3 = cnt;
                        prev_i_Start_x_addr3 = i_Start;
                        prev_i_End_x_addr3 = i;
                        prev_j_Start_x_addr3 = j_Start;
                        prev_j_End_x_addr3 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr4( i) == calAddrx_addr3(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr5( i) == calAddrx_addr3(i_Start, j_Start)) {
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
void ref_x_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr0 = -1;
    uint64_t prev_i_Start_x_addr0 = -1;
    uint64_t prev_i_End_x_addr0 = -1;
    uint64_t prev_cnt_x_addr4 = -1;
    uint64_t prev_i_Start_x_addr4 = -1;
    uint64_t prev_i_End_x_addr4 = -1;
    uint64_t prev_cnt_x_addr5 = -1;
    uint64_t prev_i_Start_x_addr5 = -1;
    uint64_t prev_i_End_x_addr5 = -1;
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

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_x_addr0 != -1) {
            if ( calAddrx_addr0( i_Start - prev_i_Start_x_addr0 + prev_i_End_x_addr0) == calAddrx_addr4(i_Start)) {
                rtHistoCal(prev_cnt_x_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr4 != -1) {
            if ( calAddrx_addr4( i_Start - prev_i_Start_x_addr4 + prev_i_End_x_addr4) == calAddrx_addr4(i_Start)) {
                rtHistoCal(prev_cnt_x_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr5 != -1) {
            if ( calAddrx_addr5( i_Start - prev_i_Start_x_addr5 + prev_i_End_x_addr5) == calAddrx_addr4(i_Start)) {
                rtHistoCal(prev_cnt_x_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr0( i) == calAddrx_addr4(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr0 = cnt;
                    prev_i_Start_x_addr0 = i_Start;
                    prev_i_End_x_addr0 = i;
                    goto EndSample;
                }
            }
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr4(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr2( i, j) == calAddrx_addr4(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr3( i, j) == calAddrx_addr4(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr4( i) == calAddrx_addr4(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr4 = cnt;
                    prev_i_Start_x_addr4 = i_Start;
                    prev_i_End_x_addr4 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr5( i) == calAddrx_addr4(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr5 = cnt;
                    prev_i_Start_x_addr5 = i_Start;
                    prev_i_End_x_addr5 = i;
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_x_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_x_addr0 = -1;
    uint64_t prev_i_Start_x_addr0 = -1;
    uint64_t prev_i_End_x_addr0 = -1;
    uint64_t prev_cnt_x_addr4 = -1;
    uint64_t prev_i_Start_x_addr4 = -1;
    uint64_t prev_i_End_x_addr4 = -1;
    uint64_t prev_cnt_x_addr5 = -1;
    uint64_t prev_i_Start_x_addr5 = -1;
    uint64_t prev_i_End_x_addr5 = -1;
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

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_x_addr0 != -1) {
            if ( calAddrx_addr0( i_Start - prev_i_Start_x_addr0 + prev_i_End_x_addr0) == calAddrx_addr5(i_Start)) {
                rtHistoCal(prev_cnt_x_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr4 != -1) {
            if ( calAddrx_addr4( i_Start - prev_i_Start_x_addr4 + prev_i_End_x_addr4) == calAddrx_addr5(i_Start)) {
                rtHistoCal(prev_cnt_x_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_x_addr5 != -1) {
            if ( calAddrx_addr5( i_Start - prev_i_Start_x_addr5 + prev_i_End_x_addr5) == calAddrx_addr5(i_Start)) {
                rtHistoCal(prev_cnt_x_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr0( i) == calAddrx_addr5(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr0 = cnt;
                    prev_i_Start_x_addr0 = i_Start;
                    prev_i_End_x_addr0 = i;
                    goto EndSample;
                }
            }
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr5(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr2( i, j) == calAddrx_addr5(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr3( i, j) == calAddrx_addr5(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr4( i) == calAddrx_addr5(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr4 = cnt;
                    prev_i_Start_x_addr4 = i_Start;
                    prev_i_End_x_addr4 = i;
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr5( i) == calAddrx_addr5(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_x_addr5 = cnt;
                    prev_i_Start_x_addr5 = i_Start;
                    prev_i_End_x_addr5 = i;
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
    ref_L_addr0();
    ref_L_addr1();
    ref_b_addr0();
    ref_x_addr0();
    ref_x_addr1();
    ref_x_addr2();
    ref_x_addr3();
    ref_x_addr4();
    ref_x_addr5();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
trisolv */ 
