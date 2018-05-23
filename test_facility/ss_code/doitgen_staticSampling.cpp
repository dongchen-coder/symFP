
 /* Start to analysis array index
Array index info
A.addr ((((r * 256) * 256) + (q * 256)) + s)
sum.addr p
C4.addr ((s * 256) + p)
sum.addr p
sum.addr p
sum.addr p
A.addr ((((r * 256) * 256) + (q * 256)) + p)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %sum
double* %A
double* %C4

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--r
--Loop Bound: (0, 256)
--Loop inc: (r + 1)
--Loop predicate: <
----q
----Loop Bound: (0, 256)
----Loop inc: (q + 1)
----Loop predicate: <
------p
------Loop Bound: (0, 256)
------Loop inc: (p + 1)
------Loop predicate: <
--------array access sum.addr p
--------s
--------Loop Bound: (0, 256)
--------Loop inc: (s + 1)
--------Loop predicate: <
----------array access A.addr ((((r * 256) * 256) + (q * 256)) + s)
----------array access C4.addr ((s * 256) + p)
----------array access sum.addr p
----------array access sum.addr p
------p
------Loop Bound: (0, 256)
------Loop inc: (p + 1)
------Loop predicate: <
--------array access sum.addr p
--------array access A.addr ((((r * 256) * 256) + (q * 256)) + p)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 2
------Sample number: 6
--------Sample number: 16
----------Sample number: 42
--------Sample number: 16
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
int calAddrsum_addr0( int r, int q, int p) {
    int result = (p) * 8 / 64;
    return result;
}
int calAddrA_addr0( int r, int q, int p, int s) {
    int result = (((((r * 256) * 256) + (q * 256)) + s)) * 8 / 64;
    return result;
}
int calAddrC4_addr0( int r, int q, int p, int s) {
    int result = (((s * 256) + p)) * 8 / 64;
    return result;
}
int calAddrsum_addr1( int r, int q, int p, int s) {
    int result = (p) * 8 / 64;
    return result;
}
int calAddrsum_addr2( int r, int q, int p, int s) {
    int result = (p) * 8 / 64;
    return result;
}
int calAddrsum_addr3( int r, int q, int p) {
    int result = (p) * 8 / 64;
    return result;
}
int calAddrA_addr1( int r, int q, int p) {
    int result = (((((r * 256) * 256) + (q * 256)) + p)) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_r_Start_A_addr0 = -1;
    uint64_t prev_r_End_A_addr0 = -1;
    uint64_t prev_q_Start_A_addr0 = -1;
    uint64_t prev_q_End_A_addr0 = -1;
    uint64_t prev_p_Start_A_addr0 = -1;
    uint64_t prev_p_End_A_addr0 = -1;
    uint64_t prev_s_Start_A_addr0 = -1;
    uint64_t prev_s_End_A_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 42;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int s_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" + std::to_string(s_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( r_Start - prev_r_Start_A_addr0 + prev_r_End_A_addr0, q_Start - prev_q_Start_A_addr0 + prev_q_End_A_addr0, p_Start - prev_p_Start_A_addr0 + prev_p_End_A_addr0, s_Start - prev_s_Start_A_addr0 + prev_s_End_A_addr0) == calAddrA_addr0(r_Start, q_Start, p_Start, s_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int rLB0 = r_Start;
        for ( int r = rLB0; r < 256; r++) {
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if (cntStart == true) cnt++;
                    {
                    int sLB3 = 0;
                    if ( r == r_Start && q == q_Start && p == p_Start ) {
                        sLB3 = s_Start;
                    }
                    for ( int s = sLB3; s < 256; s++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( r, q, p, s) == calAddrA_addr0(r_Start, q_Start, p_Start, s_Start)) {
                                rtHistoCal(cnt);
                                prev_cnt_A_addr0 = cnt;
                                prev_r_Start_A_addr0 = r_Start;
                                prev_r_End_A_addr0 = r;
                                prev_q_Start_A_addr0 = q_Start;
                                prev_q_End_A_addr0 = q;
                                prev_p_Start_A_addr0 = p_Start;
                                prev_p_End_A_addr0 = p;
                                prev_s_Start_A_addr0 = s_Start;
                                prev_s_End_A_addr0 = s;
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
                {
                int pLB4 = 0;
                for ( int p = pLB4; p < 256; p++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( r, q, p) == calAddrA_addr0(r_Start, q_Start, p_Start, s_Start)) {
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
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_r_Start_A_addr1 = -1;
    uint64_t prev_r_End_A_addr1 = -1;
    uint64_t prev_q_Start_A_addr1 = -1;
    uint64_t prev_q_End_A_addr1 = -1;
    uint64_t prev_p_Start_A_addr1 = -1;
    uint64_t prev_p_End_A_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( r_Start - prev_r_Start_A_addr1 + prev_r_End_A_addr1, q_Start - prev_q_Start_A_addr1 + prev_q_End_A_addr1, p_Start - prev_p_Start_A_addr1 + prev_p_End_A_addr1) == calAddrA_addr1(r_Start, q_Start, p_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int rLB0 = r_Start;
        for ( int r = rLB0; r < 256; r++) {
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                {
                int pLB2 = 0;
                for ( int p = pLB2; p < 256; p++) {
                    if (cntStart == true) cnt++;
                    {
                    int sLB3 = 0;
                    for ( int s = sLB3; s < 256; s++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( r, q, p, s) == calAddrA_addr1(r_Start, q_Start, p_Start)) {
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
                {
                int pLB4 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB4 = p_Start;
                }
                for ( int p = pLB4; p < 256; p++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( r, q, p) == calAddrA_addr1(r_Start, q_Start, p_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_r_Start_A_addr1 = r_Start;
                            prev_r_End_A_addr1 = r;
                            prev_q_Start_A_addr1 = q_Start;
                            prev_q_End_A_addr1 = q;
                            prev_p_Start_A_addr1 = p_Start;
                            prev_p_End_A_addr1 = p;
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
void ref_C4_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_C4_addr0 = -1;
    uint64_t prev_r_Start_C4_addr0 = -1;
    uint64_t prev_r_End_C4_addr0 = -1;
    uint64_t prev_q_Start_C4_addr0 = -1;
    uint64_t prev_q_End_C4_addr0 = -1;
    uint64_t prev_p_Start_C4_addr0 = -1;
    uint64_t prev_p_End_C4_addr0 = -1;
    uint64_t prev_s_Start_C4_addr0 = -1;
    uint64_t prev_s_End_C4_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 42;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int s_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" + std::to_string(s_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_C4_addr0 != -1) {
            if ( calAddrC4_addr0( r_Start - prev_r_Start_C4_addr0 + prev_r_End_C4_addr0, q_Start - prev_q_Start_C4_addr0 + prev_q_End_C4_addr0, p_Start - prev_p_Start_C4_addr0 + prev_p_End_C4_addr0, s_Start - prev_s_Start_C4_addr0 + prev_s_End_C4_addr0) == calAddrC4_addr0(r_Start, q_Start, p_Start, s_Start)) {
                rtHistoCal(prev_cnt_C4_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int rLB0 = r_Start;
        for ( int r = rLB0; r < 256; r++) {
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if (cntStart == true) cnt++;
                    {
                    int sLB3 = 0;
                    if ( r == r_Start && q == q_Start && p == p_Start ) {
                        sLB3 = s_Start;
                    }
                    for ( int s = sLB3; s < 256; s++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC4_addr0( r, q, p, s) == calAddrC4_addr0(r_Start, q_Start, p_Start, s_Start)) {
                                rtHistoCal(cnt);
                                prev_cnt_C4_addr0 = cnt;
                                prev_r_Start_C4_addr0 = r_Start;
                                prev_r_End_C4_addr0 = r;
                                prev_q_Start_C4_addr0 = q_Start;
                                prev_q_End_C4_addr0 = q;
                                prev_p_Start_C4_addr0 = p_Start;
                                prev_p_End_C4_addr0 = p;
                                prev_s_Start_C4_addr0 = s_Start;
                                prev_s_End_C4_addr0 = s;
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
                int pLB4 = 0;
                for ( int p = pLB4; p < 256; p++) {
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
void ref_sum_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_sum_addr0 = -1;
    uint64_t prev_r_Start_sum_addr0 = -1;
    uint64_t prev_r_End_sum_addr0 = -1;
    uint64_t prev_q_Start_sum_addr0 = -1;
    uint64_t prev_q_End_sum_addr0 = -1;
    uint64_t prev_p_Start_sum_addr0 = -1;
    uint64_t prev_p_End_sum_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_sum_addr0 != -1) {
            if ( calAddrsum_addr0( r_Start - prev_r_Start_sum_addr0 + prev_r_End_sum_addr0, q_Start - prev_q_Start_sum_addr0 + prev_q_End_sum_addr0, p_Start - prev_p_Start_sum_addr0 + prev_p_End_sum_addr0) == calAddrsum_addr0(r_Start, q_Start, p_Start)) {
                rtHistoCal(prev_cnt_sum_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int rLB0 = r_Start;
        for ( int r = rLB0; r < 256; r++) {
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr0( r, q, p) == calAddrsum_addr0(r_Start, q_Start, p_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_sum_addr0 = cnt;
                            prev_r_Start_sum_addr0 = r_Start;
                            prev_r_End_sum_addr0 = r;
                            prev_q_Start_sum_addr0 = q_Start;
                            prev_q_End_sum_addr0 = q;
                            prev_p_Start_sum_addr0 = p_Start;
                            prev_p_End_sum_addr0 = p;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    {
                    int sLB3 = 0;
                    for ( int s = sLB3; s < 256; s++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr1( r, q, p, s) == calAddrsum_addr0(r_Start, q_Start, p_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr2( r, q, p, s) == calAddrsum_addr0(r_Start, q_Start, p_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                    }
                    }
                }
                }
                {
                int pLB4 = 0;
                for ( int p = pLB4; p < 256; p++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr3( r, q, p) == calAddrsum_addr0(r_Start, q_Start, p_Start)) {
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
void ref_sum_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_sum_addr1 = -1;
    uint64_t prev_r_Start_sum_addr1 = -1;
    uint64_t prev_r_End_sum_addr1 = -1;
    uint64_t prev_q_Start_sum_addr1 = -1;
    uint64_t prev_q_End_sum_addr1 = -1;
    uint64_t prev_p_Start_sum_addr1 = -1;
    uint64_t prev_p_End_sum_addr1 = -1;
    uint64_t prev_s_Start_sum_addr1 = -1;
    uint64_t prev_s_End_sum_addr1 = -1;
    uint64_t prev_cnt_sum_addr2 = -1;
    uint64_t prev_r_Start_sum_addr2 = -1;
    uint64_t prev_r_End_sum_addr2 = -1;
    uint64_t prev_q_Start_sum_addr2 = -1;
    uint64_t prev_q_End_sum_addr2 = -1;
    uint64_t prev_p_Start_sum_addr2 = -1;
    uint64_t prev_p_End_sum_addr2 = -1;
    uint64_t prev_s_Start_sum_addr2 = -1;
    uint64_t prev_s_End_sum_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 42;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int s_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" + std::to_string(s_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_sum_addr1 != -1) {
            if ( calAddrsum_addr1( r_Start - prev_r_Start_sum_addr1 + prev_r_End_sum_addr1, q_Start - prev_q_Start_sum_addr1 + prev_q_End_sum_addr1, p_Start - prev_p_Start_sum_addr1 + prev_p_End_sum_addr1, s_Start - prev_s_Start_sum_addr1 + prev_s_End_sum_addr1) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
                rtHistoCal(prev_cnt_sum_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_sum_addr2 != -1) {
            if ( calAddrsum_addr2( r_Start - prev_r_Start_sum_addr2 + prev_r_End_sum_addr2, q_Start - prev_q_Start_sum_addr2 + prev_q_End_sum_addr2, p_Start - prev_p_Start_sum_addr2 + prev_p_End_sum_addr2, s_Start - prev_s_Start_sum_addr2 + prev_s_End_sum_addr2) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
                rtHistoCal(prev_cnt_sum_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int rLB0 = r_Start;
        for ( int r = rLB0; r < 256; r++) {
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr0( r, q, p) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    {
                    int sLB3 = 0;
                    if ( r == r_Start && q == q_Start && p == p_Start ) {
                        sLB3 = s_Start;
                    }
                    for ( int s = sLB3; s < 256; s++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr1( r, q, p, s) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
                                rtHistoCal(cnt);
                                prev_cnt_sum_addr1 = cnt;
                                prev_r_Start_sum_addr1 = r_Start;
                                prev_r_End_sum_addr1 = r;
                                prev_q_Start_sum_addr1 = q_Start;
                                prev_q_End_sum_addr1 = q;
                                prev_p_Start_sum_addr1 = p_Start;
                                prev_p_End_sum_addr1 = p;
                                prev_s_Start_sum_addr1 = s_Start;
                                prev_s_End_sum_addr1 = s;
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr2( r, q, p, s) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
                                rtHistoCal(cnt);
                                prev_cnt_sum_addr2 = cnt;
                                prev_r_Start_sum_addr2 = r_Start;
                                prev_r_End_sum_addr2 = r;
                                prev_q_Start_sum_addr2 = q_Start;
                                prev_q_End_sum_addr2 = q;
                                prev_p_Start_sum_addr2 = p_Start;
                                prev_p_End_sum_addr2 = p;
                                prev_s_Start_sum_addr2 = s_Start;
                                prev_s_End_sum_addr2 = s;
                                goto EndSample;
                            }
                        }
                    }
                    }
                }
                }
                {
                int pLB4 = 0;
                for ( int p = pLB4; p < 256; p++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr3( r, q, p) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
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
void ref_sum_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_sum_addr1 = -1;
    uint64_t prev_r_Start_sum_addr1 = -1;
    uint64_t prev_r_End_sum_addr1 = -1;
    uint64_t prev_q_Start_sum_addr1 = -1;
    uint64_t prev_q_End_sum_addr1 = -1;
    uint64_t prev_p_Start_sum_addr1 = -1;
    uint64_t prev_p_End_sum_addr1 = -1;
    uint64_t prev_s_Start_sum_addr1 = -1;
    uint64_t prev_s_End_sum_addr1 = -1;
    uint64_t prev_cnt_sum_addr2 = -1;
    uint64_t prev_r_Start_sum_addr2 = -1;
    uint64_t prev_r_End_sum_addr2 = -1;
    uint64_t prev_q_Start_sum_addr2 = -1;
    uint64_t prev_q_End_sum_addr2 = -1;
    uint64_t prev_p_Start_sum_addr2 = -1;
    uint64_t prev_p_End_sum_addr2 = -1;
    uint64_t prev_s_Start_sum_addr2 = -1;
    uint64_t prev_s_End_sum_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 42;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int s_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" + std::to_string(s_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_sum_addr1 != -1) {
            if ( calAddrsum_addr1( r_Start - prev_r_Start_sum_addr1 + prev_r_End_sum_addr1, q_Start - prev_q_Start_sum_addr1 + prev_q_End_sum_addr1, p_Start - prev_p_Start_sum_addr1 + prev_p_End_sum_addr1, s_Start - prev_s_Start_sum_addr1 + prev_s_End_sum_addr1) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
                rtHistoCal(prev_cnt_sum_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_sum_addr2 != -1) {
            if ( calAddrsum_addr2( r_Start - prev_r_Start_sum_addr2 + prev_r_End_sum_addr2, q_Start - prev_q_Start_sum_addr2 + prev_q_End_sum_addr2, p_Start - prev_p_Start_sum_addr2 + prev_p_End_sum_addr2, s_Start - prev_s_Start_sum_addr2 + prev_s_End_sum_addr2) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
                rtHistoCal(prev_cnt_sum_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int rLB0 = r_Start;
        for ( int r = rLB0; r < 256; r++) {
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr0( r, q, p) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    {
                    int sLB3 = 0;
                    if ( r == r_Start && q == q_Start && p == p_Start ) {
                        sLB3 = s_Start;
                    }
                    for ( int s = sLB3; s < 256; s++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr1( r, q, p, s) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
                                rtHistoCal(cnt);
                                prev_cnt_sum_addr1 = cnt;
                                prev_r_Start_sum_addr1 = r_Start;
                                prev_r_End_sum_addr1 = r;
                                prev_q_Start_sum_addr1 = q_Start;
                                prev_q_End_sum_addr1 = q;
                                prev_p_Start_sum_addr1 = p_Start;
                                prev_p_End_sum_addr1 = p;
                                prev_s_Start_sum_addr1 = s_Start;
                                prev_s_End_sum_addr1 = s;
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr2( r, q, p, s) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
                                rtHistoCal(cnt);
                                prev_cnt_sum_addr2 = cnt;
                                prev_r_Start_sum_addr2 = r_Start;
                                prev_r_End_sum_addr2 = r;
                                prev_q_Start_sum_addr2 = q_Start;
                                prev_q_End_sum_addr2 = q;
                                prev_p_Start_sum_addr2 = p_Start;
                                prev_p_End_sum_addr2 = p;
                                prev_s_Start_sum_addr2 = s_Start;
                                prev_s_End_sum_addr2 = s;
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                    }
                    }
                }
                }
                {
                int pLB4 = 0;
                for ( int p = pLB4; p < 256; p++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr3( r, q, p) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
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
void ref_sum_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_sum_addr3 = -1;
    uint64_t prev_r_Start_sum_addr3 = -1;
    uint64_t prev_r_End_sum_addr3 = -1;
    uint64_t prev_q_Start_sum_addr3 = -1;
    uint64_t prev_q_End_sum_addr3 = -1;
    uint64_t prev_p_Start_sum_addr3 = -1;
    uint64_t prev_p_End_sum_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_sum_addr3 != -1) {
            if ( calAddrsum_addr3( r_Start - prev_r_Start_sum_addr3 + prev_r_End_sum_addr3, q_Start - prev_q_Start_sum_addr3 + prev_q_End_sum_addr3, p_Start - prev_p_Start_sum_addr3 + prev_p_End_sum_addr3) == calAddrsum_addr3(r_Start, q_Start, p_Start)) {
                rtHistoCal(prev_cnt_sum_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int rLB0 = r_Start;
        for ( int r = rLB0; r < 256; r++) {
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                {
                int pLB2 = 0;
                for ( int p = pLB2; p < 256; p++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr0( r, q, p) == calAddrsum_addr3(r_Start, q_Start, p_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    {
                    int sLB3 = 0;
                    for ( int s = sLB3; s < 256; s++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr1( r, q, p, s) == calAddrsum_addr3(r_Start, q_Start, p_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr2( r, q, p, s) == calAddrsum_addr3(r_Start, q_Start, p_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                    }
                    }
                }
                }
                {
                int pLB4 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB4 = p_Start;
                }
                for ( int p = pLB4; p < 256; p++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr3( r, q, p) == calAddrsum_addr3(r_Start, q_Start, p_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_sum_addr3 = cnt;
                            prev_r_Start_sum_addr3 = r_Start;
                            prev_r_End_sum_addr3 = r;
                            prev_q_Start_sum_addr3 = q_Start;
                            prev_q_End_sum_addr3 = q;
                            prev_p_Start_sum_addr3 = p_Start;
                            prev_p_End_sum_addr3 = p;
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
int main() {
    ref_A_addr0();
    ref_A_addr1();
    ref_C4_addr0();
    ref_sum_addr0();
    ref_sum_addr1();
    ref_sum_addr2();
    ref_sum_addr3();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
doitgen */ 
