
 /* Start to analysis array index
Array index info
q.addr i
s.addr i
s.addr j
r.addr i
A.addr ((i * 1024) + j)
s.addr j
q.addr i
A.addr ((i * 1024) + j)
p.addr j
q.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %nx
i32 %ny
double* %A
double* %r
double* %s
double* %p
double* %q

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----array access s.addr i
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----array access q.addr i
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access s.addr j
------array access r.addr i
------array access A.addr ((i * 1024) + j)
------array access s.addr j
------array access q.addr i
------array access A.addr ((i * 1024) + j)
------array access p.addr j
------array access q.addr i

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 1
----Sample number: 1
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
int calAddrs_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrq_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrs_addr1( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrr_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrA_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrs_addr2( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrq_addr1( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrA_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrp_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrq_addr2( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
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
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0) == calAddrA_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1) == calAddrA_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr0 = cnt;
                        prev_i_Start_A_addr0 = i_Start;
                        prev_i_End_A_addr0 = i;
                        prev_j_Start_A_addr0 = j_Start;
                        prev_j_End_A_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr1 = cnt;
                        prev_i_Start_A_addr1 = i_Start;
                        prev_i_End_A_addr1 = i;
                        prev_j_Start_A_addr1 = j_Start;
                        prev_j_End_A_addr1 = j;
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
void ref_A_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
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
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0) == calAddrA_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1) == calAddrA_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( i, j) == calAddrA_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr0 = cnt;
                        prev_i_Start_A_addr0 = i_Start;
                        prev_i_End_A_addr0 = i;
                        prev_j_Start_A_addr0 = j_Start;
                        prev_j_End_A_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( i, j) == calAddrA_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr1 = cnt;
                        prev_i_Start_A_addr1 = i_Start;
                        prev_i_End_A_addr1 = i;
                        prev_j_Start_A_addr1 = j_Start;
                        prev_j_End_A_addr1 = j;
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
void ref_p_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr0 = -1;
    uint64_t prev_i_Start_p_addr0 = -1;
    uint64_t prev_i_End_p_addr0 = -1;
    uint64_t prev_j_Start_p_addr0 = -1;
    uint64_t prev_j_End_p_addr0 = -1;
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
        if ( prev_cnt_p_addr0 != -1) {
            if ( calAddrp_addr0( i_Start - prev_i_Start_p_addr0 + prev_i_End_p_addr0, j_Start - prev_j_Start_p_addr0 + prev_j_End_p_addr0) == calAddrp_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( i, j) == calAddrp_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_p_addr0 = cnt;
                        prev_i_Start_p_addr0 = i_Start;
                        prev_i_End_p_addr0 = i;
                        prev_j_Start_p_addr0 = j_Start;
                        prev_j_End_p_addr0 = j;
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
void ref_q_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr0 = -1;
    uint64_t prev_i_Start_q_addr0 = -1;
    uint64_t prev_i_End_q_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr0 != -1) {
            if ( calAddrq_addr0( i_Start - prev_i_Start_q_addr0 + prev_i_End_q_addr0) == calAddrq_addr0(i_Start)) {
                rtHistoCal(prev_cnt_q_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr0( i) == calAddrq_addr0(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_q_addr0 = cnt;
                    prev_i_Start_q_addr0 = i_Start;
                    prev_i_End_q_addr0 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr1( i, j) == calAddrq_addr0(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr2( i, j) == calAddrq_addr0(i_Start)) {
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
void ref_q_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr1 = -1;
    uint64_t prev_i_Start_q_addr1 = -1;
    uint64_t prev_i_End_q_addr1 = -1;
    uint64_t prev_j_Start_q_addr1 = -1;
    uint64_t prev_j_End_q_addr1 = -1;
    uint64_t prev_cnt_q_addr2 = -1;
    uint64_t prev_i_Start_q_addr2 = -1;
    uint64_t prev_i_End_q_addr2 = -1;
    uint64_t prev_j_Start_q_addr2 = -1;
    uint64_t prev_j_End_q_addr2 = -1;
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
        if ( prev_cnt_q_addr1 != -1) {
            if ( calAddrq_addr1( i_Start - prev_i_Start_q_addr1 + prev_i_End_q_addr1, j_Start - prev_j_Start_q_addr1 + prev_j_End_q_addr1) == calAddrq_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr2 != -1) {
            if ( calAddrq_addr2( i_Start - prev_i_Start_q_addr2 + prev_i_End_q_addr2, j_Start - prev_j_Start_q_addr2 + prev_j_End_q_addr2) == calAddrq_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr0( i) == calAddrq_addr1(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr1( i, j) == calAddrq_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_q_addr1 = cnt;
                        prev_i_Start_q_addr1 = i_Start;
                        prev_i_End_q_addr1 = i;
                        prev_j_Start_q_addr1 = j_Start;
                        prev_j_End_q_addr1 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr2( i, j) == calAddrq_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_q_addr2 = cnt;
                        prev_i_Start_q_addr2 = i_Start;
                        prev_i_End_q_addr2 = i;
                        prev_j_Start_q_addr2 = j_Start;
                        prev_j_End_q_addr2 = j;
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
void ref_q_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr1 = -1;
    uint64_t prev_i_Start_q_addr1 = -1;
    uint64_t prev_i_End_q_addr1 = -1;
    uint64_t prev_j_Start_q_addr1 = -1;
    uint64_t prev_j_End_q_addr1 = -1;
    uint64_t prev_cnt_q_addr2 = -1;
    uint64_t prev_i_Start_q_addr2 = -1;
    uint64_t prev_i_End_q_addr2 = -1;
    uint64_t prev_j_Start_q_addr2 = -1;
    uint64_t prev_j_End_q_addr2 = -1;
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
        if ( prev_cnt_q_addr1 != -1) {
            if ( calAddrq_addr1( i_Start - prev_i_Start_q_addr1 + prev_i_End_q_addr1, j_Start - prev_j_Start_q_addr1 + prev_j_End_q_addr1) == calAddrq_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr2 != -1) {
            if ( calAddrq_addr2( i_Start - prev_i_Start_q_addr2 + prev_i_End_q_addr2, j_Start - prev_j_Start_q_addr2 + prev_j_End_q_addr2) == calAddrq_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr0( i) == calAddrq_addr2(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr1( i, j) == calAddrq_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_q_addr1 = cnt;
                        prev_i_Start_q_addr1 = i_Start;
                        prev_i_End_q_addr1 = i;
                        prev_j_Start_q_addr1 = j_Start;
                        prev_j_End_q_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr2( i, j) == calAddrq_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_q_addr2 = cnt;
                        prev_i_Start_q_addr2 = i_Start;
                        prev_i_End_q_addr2 = i;
                        prev_j_Start_q_addr2 = j_Start;
                        prev_j_End_q_addr2 = j;
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
void ref_r_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_r_addr0 = -1;
    uint64_t prev_i_Start_r_addr0 = -1;
    uint64_t prev_i_End_r_addr0 = -1;
    uint64_t prev_j_Start_r_addr0 = -1;
    uint64_t prev_j_End_r_addr0 = -1;
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
        if ( prev_cnt_r_addr0 != -1) {
            if ( calAddrr_addr0( i_Start - prev_i_Start_r_addr0 + prev_i_End_r_addr0, j_Start - prev_j_Start_r_addr0 + prev_j_End_r_addr0) == calAddrr_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_r_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrr_addr0( i, j) == calAddrr_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_r_addr0 = cnt;
                        prev_i_Start_r_addr0 = i_Start;
                        prev_i_End_r_addr0 = i;
                        prev_j_Start_r_addr0 = j_Start;
                        prev_j_End_r_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
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
EndSample:
        s++;
        }
}
void ref_s_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_s_addr0 = -1;
    uint64_t prev_i_Start_s_addr0 = -1;
    uint64_t prev_i_End_s_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_s_addr0 != -1) {
            if ( calAddrs_addr0( i_Start - prev_i_Start_s_addr0 + prev_i_End_s_addr0) == calAddrs_addr0(i_Start)) {
                rtHistoCal(prev_cnt_s_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrs_addr0( i) == calAddrs_addr0(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_s_addr0 = cnt;
                    prev_i_Start_s_addr0 = i_Start;
                    prev_i_End_s_addr0 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
        }
        }
        {
        int iLB1 = 0;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr1( i, j) == calAddrs_addr0(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr2( i, j) == calAddrs_addr0(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
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
void ref_s_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_s_addr1 = -1;
    uint64_t prev_i_Start_s_addr1 = -1;
    uint64_t prev_i_End_s_addr1 = -1;
    uint64_t prev_j_Start_s_addr1 = -1;
    uint64_t prev_j_End_s_addr1 = -1;
    uint64_t prev_cnt_s_addr2 = -1;
    uint64_t prev_i_Start_s_addr2 = -1;
    uint64_t prev_i_End_s_addr2 = -1;
    uint64_t prev_j_Start_s_addr2 = -1;
    uint64_t prev_j_End_s_addr2 = -1;
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
        if ( prev_cnt_s_addr1 != -1) {
            if ( calAddrs_addr1( i_Start - prev_i_Start_s_addr1 + prev_i_End_s_addr1, j_Start - prev_j_Start_s_addr1 + prev_j_End_s_addr1) == calAddrs_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_s_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_s_addr2 != -1) {
            if ( calAddrs_addr2( i_Start - prev_i_Start_s_addr2 + prev_i_End_s_addr2, j_Start - prev_j_Start_s_addr2 + prev_j_End_s_addr2) == calAddrs_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_s_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr1( i, j) == calAddrs_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_s_addr1 = cnt;
                        prev_i_Start_s_addr1 = i_Start;
                        prev_i_End_s_addr1 = i;
                        prev_j_Start_s_addr1 = j_Start;
                        prev_j_End_s_addr1 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr2( i, j) == calAddrs_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_s_addr2 = cnt;
                        prev_i_Start_s_addr2 = i_Start;
                        prev_i_End_s_addr2 = i;
                        prev_j_Start_s_addr2 = j_Start;
                        prev_j_End_s_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
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
void ref_s_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_s_addr1 = -1;
    uint64_t prev_i_Start_s_addr1 = -1;
    uint64_t prev_i_End_s_addr1 = -1;
    uint64_t prev_j_Start_s_addr1 = -1;
    uint64_t prev_j_End_s_addr1 = -1;
    uint64_t prev_cnt_s_addr2 = -1;
    uint64_t prev_i_Start_s_addr2 = -1;
    uint64_t prev_i_End_s_addr2 = -1;
    uint64_t prev_j_Start_s_addr2 = -1;
    uint64_t prev_j_End_s_addr2 = -1;
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
        if ( prev_cnt_s_addr1 != -1) {
            if ( calAddrs_addr1( i_Start - prev_i_Start_s_addr1 + prev_i_End_s_addr1, j_Start - prev_j_Start_s_addr1 + prev_j_End_s_addr1) == calAddrs_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_s_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_s_addr2 != -1) {
            if ( calAddrs_addr2( i_Start - prev_i_Start_s_addr2 + prev_i_End_s_addr2, j_Start - prev_j_Start_s_addr2 + prev_j_End_s_addr2) == calAddrs_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_s_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr1( i, j) == calAddrs_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_s_addr1 = cnt;
                        prev_i_Start_s_addr1 = i_Start;
                        prev_i_End_s_addr1 = i;
                        prev_j_Start_s_addr1 = j_Start;
                        prev_j_End_s_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr2( i, j) == calAddrs_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_s_addr2 = cnt;
                        prev_i_Start_s_addr2 = i_Start;
                        prev_i_End_s_addr2 = i;
                        prev_j_Start_s_addr2 = j_Start;
                        prev_j_End_s_addr2 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
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
int main() {
    ref_A_addr0();
    ref_A_addr1();
    ref_p_addr0();
    ref_q_addr0();
    ref_q_addr1();
    ref_q_addr2();
    ref_r_addr0();
    ref_s_addr0();
    ref_s_addr1();
    ref_s_addr2();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: bicg_cpu */ 
