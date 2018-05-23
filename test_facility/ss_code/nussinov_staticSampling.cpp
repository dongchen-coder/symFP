
 /* Start to analysis array index
Array index info
table.addr ((i * 1024) + j)
table.addr (((i * 1024) + j) - 1)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr (((i + 1) * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((((i + 1) * 1024) + j) - 1)
seq.addr i
seq.addr j
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((((i + 1) * 1024) + j) - 1)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + k)
table.addr (((k + 1) * 1024) + j)
table.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %table
double* %seq

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (1023, 0)
--Loop inc: (i + -1)
--Loop predicate: >=
----j
----Loop Bound: ((i + 1), 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access table.addr ((i * 1024) + j)
------array access table.addr (((i * 1024) + j) - 1)
------array access table.addr ((i * 1024) + j)
------array access table.addr ((i * 1024) + j)
------array access table.addr (((i + 1) * 1024) + j)
------array access table.addr ((i * 1024) + j)
------array access table.addr ((i * 1024) + j)
------array access table.addr ((((i + 1) * 1024) + j) - 1)
------array access seq.addr i
------array access seq.addr j
------array access table.addr ((i * 1024) + j)
------array access table.addr ((i * 1024) + j)
------array access table.addr ((((i + 1) * 1024) + j) - 1)
------array access table.addr ((i * 1024) + j)
------k
------Loop Bound: ((i + 1), j)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access table.addr ((i * 1024) + j)
--------array access table.addr ((i * 1024) + k)
--------array access table.addr (((k + 1) * 1024) + j)
--------array access table.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 1023 1024 
Dump stride: -1 1 
init counter: 1023 1024 1024 
Dump stride: -1 1 1 
Dump tree:
----Sample number: 10
------Sample number: 52
--------Sample number: 178
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
int calAddrtable_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr1( int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrtable_addr2( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr3( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr4( int i, int j) {
    int result = ((((i + 1) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr6( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr7( int i, int j) {
    int result = (((((i + 1) * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrseq_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrseq_addr1( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
int calAddrtable_addr8( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr9( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr10( int i, int j) {
    int result = (((((i + 1) * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrtable_addr11( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr12( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr13( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrtable_addr14( int i, int j, int k) {
    int result = ((((k + 1) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrtable_addr15( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_seq_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_seq_addr0 = -1;
    uint64_t prev_i_Start_seq_addr0 = -1;
    uint64_t prev_i_End_seq_addr0 = -1;
    uint64_t prev_j_Start_seq_addr0 = -1;
    uint64_t prev_j_End_seq_addr0 = -1;
    uint64_t prev_cnt_seq_addr1 = -1;
    uint64_t prev_i_Start_seq_addr1 = -1;
    uint64_t prev_i_End_seq_addr1 = -1;
    uint64_t prev_j_Start_seq_addr1 = -1;
    uint64_t prev_j_End_seq_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_seq_addr0 != -1) {
            if ( calAddrseq_addr0( i_Start - prev_i_Start_seq_addr0 + prev_i_End_seq_addr0, j_Start - prev_j_Start_seq_addr0 + prev_j_End_seq_addr0) == calAddrseq_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_seq_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_seq_addr1 != -1) {
            if ( calAddrseq_addr1( i_Start - prev_i_Start_seq_addr1 + prev_i_End_seq_addr1, j_Start - prev_j_Start_seq_addr1 + prev_j_End_seq_addr1) == calAddrseq_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_seq_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrseq_addr0( i, j) == calAddrseq_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_seq_addr0 = cnt;
                        prev_i_Start_seq_addr0 = i_Start;
                        prev_i_End_seq_addr0 = i;
                        prev_j_Start_seq_addr0 = j_Start;
                        prev_j_End_seq_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrseq_addr1( i, j) == calAddrseq_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_seq_addr1 = cnt;
                        prev_i_Start_seq_addr1 = i_Start;
                        prev_i_End_seq_addr1 = i;
                        prev_j_Start_seq_addr1 = j_Start;
                        prev_j_End_seq_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
void ref_seq_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_seq_addr0 = -1;
    uint64_t prev_i_Start_seq_addr0 = -1;
    uint64_t prev_i_End_seq_addr0 = -1;
    uint64_t prev_j_Start_seq_addr0 = -1;
    uint64_t prev_j_End_seq_addr0 = -1;
    uint64_t prev_cnt_seq_addr1 = -1;
    uint64_t prev_i_Start_seq_addr1 = -1;
    uint64_t prev_i_End_seq_addr1 = -1;
    uint64_t prev_j_Start_seq_addr1 = -1;
    uint64_t prev_j_End_seq_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_seq_addr0 != -1) {
            if ( calAddrseq_addr0( i_Start - prev_i_Start_seq_addr0 + prev_i_End_seq_addr0, j_Start - prev_j_Start_seq_addr0 + prev_j_End_seq_addr0) == calAddrseq_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_seq_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_seq_addr1 != -1) {
            if ( calAddrseq_addr1( i_Start - prev_i_Start_seq_addr1 + prev_i_End_seq_addr1, j_Start - prev_j_Start_seq_addr1 + prev_j_End_seq_addr1) == calAddrseq_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_seq_addr1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrseq_addr0( i, j) == calAddrseq_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_seq_addr0 = cnt;
                        prev_i_Start_seq_addr0 = i_Start;
                        prev_i_End_seq_addr0 = i;
                        prev_j_Start_seq_addr0 = j_Start;
                        prev_j_End_seq_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrseq_addr1( i, j) == calAddrseq_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_seq_addr1 = cnt;
                        prev_i_Start_seq_addr1 = i_Start;
                        prev_i_End_seq_addr1 = i;
                        prev_j_Start_seq_addr1 = j_Start;
                        prev_j_End_seq_addr1 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
void ref_table_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr0(i_Start, j_Start)) {
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
void ref_table_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr1(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr1(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr1(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr1(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr1(i_Start, j_Start)) {
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
void ref_table_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr2(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr2(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr2(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr2(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr2(i_Start, j_Start)) {
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
void ref_table_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr3(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr3(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr3(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr3(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr3(i_Start, j_Start)) {
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
void ref_table_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr4(i_Start, j_Start)) {
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
void ref_table_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr5(i_Start, j_Start)) {
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
void ref_table_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr6(i_Start, j_Start)) {
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
void ref_table_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr7(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr7(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr7(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr7(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr7(i_Start, j_Start)) {
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
void ref_table_addr8() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr8(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr8(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr8(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr8(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr8(i_Start, j_Start)) {
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
void ref_table_addr9() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr9(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr9(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr9(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr9(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr9(i_Start, j_Start)) {
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
void ref_table_addr10() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr10(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr10(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr10(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr10(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr10(i_Start, j_Start)) {
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
void ref_table_addr11() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr0 = -1;
    uint64_t prev_i_Start_table_addr0 = -1;
    uint64_t prev_i_End_table_addr0 = -1;
    uint64_t prev_j_Start_table_addr0 = -1;
    uint64_t prev_j_End_table_addr0 = -1;
    uint64_t prev_cnt_table_addr1 = -1;
    uint64_t prev_i_Start_table_addr1 = -1;
    uint64_t prev_i_End_table_addr1 = -1;
    uint64_t prev_j_Start_table_addr1 = -1;
    uint64_t prev_j_End_table_addr1 = -1;
    uint64_t prev_cnt_table_addr2 = -1;
    uint64_t prev_i_Start_table_addr2 = -1;
    uint64_t prev_i_End_table_addr2 = -1;
    uint64_t prev_j_Start_table_addr2 = -1;
    uint64_t prev_j_End_table_addr2 = -1;
    uint64_t prev_cnt_table_addr3 = -1;
    uint64_t prev_i_Start_table_addr3 = -1;
    uint64_t prev_i_End_table_addr3 = -1;
    uint64_t prev_j_Start_table_addr3 = -1;
    uint64_t prev_j_End_table_addr3 = -1;
    uint64_t prev_cnt_table_addr4 = -1;
    uint64_t prev_i_Start_table_addr4 = -1;
    uint64_t prev_i_End_table_addr4 = -1;
    uint64_t prev_j_Start_table_addr4 = -1;
    uint64_t prev_j_End_table_addr4 = -1;
    uint64_t prev_cnt_table_addr5 = -1;
    uint64_t prev_i_Start_table_addr5 = -1;
    uint64_t prev_i_End_table_addr5 = -1;
    uint64_t prev_j_Start_table_addr5 = -1;
    uint64_t prev_j_End_table_addr5 = -1;
    uint64_t prev_cnt_table_addr6 = -1;
    uint64_t prev_i_Start_table_addr6 = -1;
    uint64_t prev_i_End_table_addr6 = -1;
    uint64_t prev_j_Start_table_addr6 = -1;
    uint64_t prev_j_End_table_addr6 = -1;
    uint64_t prev_cnt_table_addr7 = -1;
    uint64_t prev_i_Start_table_addr7 = -1;
    uint64_t prev_i_End_table_addr7 = -1;
    uint64_t prev_j_Start_table_addr7 = -1;
    uint64_t prev_j_End_table_addr7 = -1;
    uint64_t prev_cnt_table_addr8 = -1;
    uint64_t prev_i_Start_table_addr8 = -1;
    uint64_t prev_i_End_table_addr8 = -1;
    uint64_t prev_j_Start_table_addr8 = -1;
    uint64_t prev_j_End_table_addr8 = -1;
    uint64_t prev_cnt_table_addr9 = -1;
    uint64_t prev_i_Start_table_addr9 = -1;
    uint64_t prev_i_End_table_addr9 = -1;
    uint64_t prev_j_Start_table_addr9 = -1;
    uint64_t prev_j_End_table_addr9 = -1;
    uint64_t prev_cnt_table_addr10 = -1;
    uint64_t prev_i_Start_table_addr10 = -1;
    uint64_t prev_i_End_table_addr10 = -1;
    uint64_t prev_j_Start_table_addr10 = -1;
    uint64_t prev_j_End_table_addr10 = -1;
    uint64_t prev_cnt_table_addr11 = -1;
    uint64_t prev_i_Start_table_addr11 = -1;
    uint64_t prev_i_End_table_addr11 = -1;
    uint64_t prev_j_Start_table_addr11 = -1;
    uint64_t prev_j_End_table_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr0 != -1) {
            if ( calAddrtable_addr0( i_Start - prev_i_Start_table_addr0 + prev_i_End_table_addr0, j_Start - prev_j_Start_table_addr0 + prev_j_End_table_addr0) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr1 != -1) {
            if ( calAddrtable_addr1( i_Start - prev_i_Start_table_addr1 + prev_i_End_table_addr1, j_Start - prev_j_Start_table_addr1 + prev_j_End_table_addr1) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr2 != -1) {
            if ( calAddrtable_addr2( i_Start - prev_i_Start_table_addr2 + prev_i_End_table_addr2, j_Start - prev_j_Start_table_addr2 + prev_j_End_table_addr2) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr3 != -1) {
            if ( calAddrtable_addr3( i_Start - prev_i_Start_table_addr3 + prev_i_End_table_addr3, j_Start - prev_j_Start_table_addr3 + prev_j_End_table_addr3) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr4 != -1) {
            if ( calAddrtable_addr4( i_Start - prev_i_Start_table_addr4 + prev_i_End_table_addr4, j_Start - prev_j_Start_table_addr4 + prev_j_End_table_addr4) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr5 != -1) {
            if ( calAddrtable_addr5( i_Start - prev_i_Start_table_addr5 + prev_i_End_table_addr5, j_Start - prev_j_Start_table_addr5 + prev_j_End_table_addr5) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr6 != -1) {
            if ( calAddrtable_addr6( i_Start - prev_i_Start_table_addr6 + prev_i_End_table_addr6, j_Start - prev_j_Start_table_addr6 + prev_j_End_table_addr6) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr7 != -1) {
            if ( calAddrtable_addr7( i_Start - prev_i_Start_table_addr7 + prev_i_End_table_addr7, j_Start - prev_j_Start_table_addr7 + prev_j_End_table_addr7) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr8 != -1) {
            if ( calAddrtable_addr8( i_Start - prev_i_Start_table_addr8 + prev_i_End_table_addr8, j_Start - prev_j_Start_table_addr8 + prev_j_End_table_addr8) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr9 != -1) {
            if ( calAddrtable_addr9( i_Start - prev_i_Start_table_addr9 + prev_i_End_table_addr9, j_Start - prev_j_Start_table_addr9 + prev_j_End_table_addr9) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr10 != -1) {
            if ( calAddrtable_addr10( i_Start - prev_i_Start_table_addr10 + prev_i_End_table_addr10, j_Start - prev_j_Start_table_addr10 + prev_j_End_table_addr10) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr10);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr11 != -1) {
            if ( calAddrtable_addr11( i_Start - prev_i_Start_table_addr11 + prev_i_End_table_addr11, j_Start - prev_j_Start_table_addr11 + prev_j_End_table_addr11) == calAddrtable_addr11(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_table_addr11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr0 = cnt;
                        prev_i_Start_table_addr0 = i_Start;
                        prev_i_End_table_addr0 = i;
                        prev_j_Start_table_addr0 = j_Start;
                        prev_j_End_table_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr1 = cnt;
                        prev_i_Start_table_addr1 = i_Start;
                        prev_i_End_table_addr1 = i;
                        prev_j_Start_table_addr1 = j_Start;
                        prev_j_End_table_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr2 = cnt;
                        prev_i_Start_table_addr2 = i_Start;
                        prev_i_End_table_addr2 = i;
                        prev_j_Start_table_addr2 = j_Start;
                        prev_j_End_table_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr3 = cnt;
                        prev_i_Start_table_addr3 = i_Start;
                        prev_i_End_table_addr3 = i;
                        prev_j_Start_table_addr3 = j_Start;
                        prev_j_End_table_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr4 = cnt;
                        prev_i_Start_table_addr4 = i_Start;
                        prev_i_End_table_addr4 = i;
                        prev_j_Start_table_addr4 = j_Start;
                        prev_j_End_table_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr5 = cnt;
                        prev_i_Start_table_addr5 = i_Start;
                        prev_i_End_table_addr5 = i;
                        prev_j_Start_table_addr5 = j_Start;
                        prev_j_End_table_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr6 = cnt;
                        prev_i_Start_table_addr6 = i_Start;
                        prev_i_End_table_addr6 = i;
                        prev_j_Start_table_addr6 = j_Start;
                        prev_j_End_table_addr6 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr7 = cnt;
                        prev_i_Start_table_addr7 = i_Start;
                        prev_i_End_table_addr7 = i;
                        prev_j_Start_table_addr7 = j_Start;
                        prev_j_End_table_addr7 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr8 = cnt;
                        prev_i_Start_table_addr8 = i_Start;
                        prev_i_End_table_addr8 = i;
                        prev_j_Start_table_addr8 = j_Start;
                        prev_j_End_table_addr8 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr9 = cnt;
                        prev_i_Start_table_addr9 = i_Start;
                        prev_i_End_table_addr9 = i;
                        prev_j_Start_table_addr9 = j_Start;
                        prev_j_End_table_addr9 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr10 = cnt;
                        prev_i_Start_table_addr10 = i_Start;
                        prev_i_End_table_addr10 = i;
                        prev_j_Start_table_addr10 = j_Start;
                        prev_j_End_table_addr10 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_table_addr11 = cnt;
                        prev_i_Start_table_addr11 = i_Start;
                        prev_i_End_table_addr11 = i;
                        prev_j_Start_table_addr11 = j_Start;
                        prev_j_End_table_addr11 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr11(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr11(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr11(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr11(i_Start, j_Start)) {
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
void ref_table_addr12() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr12 = -1;
    uint64_t prev_i_Start_table_addr12 = -1;
    uint64_t prev_i_End_table_addr12 = -1;
    uint64_t prev_j_Start_table_addr12 = -1;
    uint64_t prev_j_End_table_addr12 = -1;
    uint64_t prev_k_Start_table_addr12 = -1;
    uint64_t prev_k_End_table_addr12 = -1;
    uint64_t prev_cnt_table_addr13 = -1;
    uint64_t prev_i_Start_table_addr13 = -1;
    uint64_t prev_i_End_table_addr13 = -1;
    uint64_t prev_j_Start_table_addr13 = -1;
    uint64_t prev_j_End_table_addr13 = -1;
    uint64_t prev_k_Start_table_addr13 = -1;
    uint64_t prev_k_End_table_addr13 = -1;
    uint64_t prev_cnt_table_addr14 = -1;
    uint64_t prev_i_Start_table_addr14 = -1;
    uint64_t prev_i_End_table_addr14 = -1;
    uint64_t prev_j_Start_table_addr14 = -1;
    uint64_t prev_j_End_table_addr14 = -1;
    uint64_t prev_k_Start_table_addr14 = -1;
    uint64_t prev_k_End_table_addr14 = -1;
    uint64_t prev_cnt_table_addr15 = -1;
    uint64_t prev_i_Start_table_addr15 = -1;
    uint64_t prev_i_End_table_addr15 = -1;
    uint64_t prev_j_Start_table_addr15 = -1;
    uint64_t prev_j_End_table_addr15 = -1;
    uint64_t prev_k_Start_table_addr15 = -1;
    uint64_t prev_k_End_table_addr15 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        if ( (j_Start - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr12 != -1) {
            if ( calAddrtable_addr12( i_Start - prev_i_Start_table_addr12 + prev_i_End_table_addr12, j_Start - prev_j_Start_table_addr12 + prev_j_End_table_addr12, k_Start - prev_k_Start_table_addr12 + prev_k_End_table_addr12) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr12);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr13 != -1) {
            if ( calAddrtable_addr13( i_Start - prev_i_Start_table_addr13 + prev_i_End_table_addr13, j_Start - prev_j_Start_table_addr13 + prev_j_End_table_addr13, k_Start - prev_k_Start_table_addr13 + prev_k_End_table_addr13) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr13);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr14 != -1) {
            if ( calAddrtable_addr14( i_Start - prev_i_Start_table_addr14 + prev_i_End_table_addr14, j_Start - prev_j_Start_table_addr14 + prev_j_End_table_addr14, k_Start - prev_k_Start_table_addr14 + prev_k_End_table_addr14) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr14);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr15 != -1) {
            if ( calAddrtable_addr15( i_Start - prev_i_Start_table_addr15 + prev_i_End_table_addr15, j_Start - prev_j_Start_table_addr15 + prev_j_End_table_addr15, k_Start - prev_k_Start_table_addr15 + prev_k_End_table_addr15) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr15);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr12 = cnt;
                            prev_i_Start_table_addr12 = i_Start;
                            prev_i_End_table_addr12 = i;
                            prev_j_Start_table_addr12 = j_Start;
                            prev_j_End_table_addr12 = j;
                            prev_k_Start_table_addr12 = k_Start;
                            prev_k_End_table_addr12 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr13 = cnt;
                            prev_i_Start_table_addr13 = i_Start;
                            prev_i_End_table_addr13 = i;
                            prev_j_Start_table_addr13 = j_Start;
                            prev_j_End_table_addr13 = j;
                            prev_k_Start_table_addr13 = k_Start;
                            prev_k_End_table_addr13 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr14 = cnt;
                            prev_i_Start_table_addr14 = i_Start;
                            prev_i_End_table_addr14 = i;
                            prev_j_Start_table_addr14 = j_Start;
                            prev_j_End_table_addr14 = j;
                            prev_k_Start_table_addr14 = k_Start;
                            prev_k_End_table_addr14 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr12(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr15 = cnt;
                            prev_i_Start_table_addr15 = i_Start;
                            prev_i_End_table_addr15 = i;
                            prev_j_Start_table_addr15 = j_Start;
                            prev_j_End_table_addr15 = j;
                            prev_k_Start_table_addr15 = k_Start;
                            prev_k_End_table_addr15 = k;
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
void ref_table_addr13() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr12 = -1;
    uint64_t prev_i_Start_table_addr12 = -1;
    uint64_t prev_i_End_table_addr12 = -1;
    uint64_t prev_j_Start_table_addr12 = -1;
    uint64_t prev_j_End_table_addr12 = -1;
    uint64_t prev_k_Start_table_addr12 = -1;
    uint64_t prev_k_End_table_addr12 = -1;
    uint64_t prev_cnt_table_addr13 = -1;
    uint64_t prev_i_Start_table_addr13 = -1;
    uint64_t prev_i_End_table_addr13 = -1;
    uint64_t prev_j_Start_table_addr13 = -1;
    uint64_t prev_j_End_table_addr13 = -1;
    uint64_t prev_k_Start_table_addr13 = -1;
    uint64_t prev_k_End_table_addr13 = -1;
    uint64_t prev_cnt_table_addr14 = -1;
    uint64_t prev_i_Start_table_addr14 = -1;
    uint64_t prev_i_End_table_addr14 = -1;
    uint64_t prev_j_Start_table_addr14 = -1;
    uint64_t prev_j_End_table_addr14 = -1;
    uint64_t prev_k_Start_table_addr14 = -1;
    uint64_t prev_k_End_table_addr14 = -1;
    uint64_t prev_cnt_table_addr15 = -1;
    uint64_t prev_i_Start_table_addr15 = -1;
    uint64_t prev_i_End_table_addr15 = -1;
    uint64_t prev_j_Start_table_addr15 = -1;
    uint64_t prev_j_End_table_addr15 = -1;
    uint64_t prev_k_Start_table_addr15 = -1;
    uint64_t prev_k_End_table_addr15 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        if ( (j_Start - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr12 != -1) {
            if ( calAddrtable_addr12( i_Start - prev_i_Start_table_addr12 + prev_i_End_table_addr12, j_Start - prev_j_Start_table_addr12 + prev_j_End_table_addr12, k_Start - prev_k_Start_table_addr12 + prev_k_End_table_addr12) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr12);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr13 != -1) {
            if ( calAddrtable_addr13( i_Start - prev_i_Start_table_addr13 + prev_i_End_table_addr13, j_Start - prev_j_Start_table_addr13 + prev_j_End_table_addr13, k_Start - prev_k_Start_table_addr13 + prev_k_End_table_addr13) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr13);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr14 != -1) {
            if ( calAddrtable_addr14( i_Start - prev_i_Start_table_addr14 + prev_i_End_table_addr14, j_Start - prev_j_Start_table_addr14 + prev_j_End_table_addr14, k_Start - prev_k_Start_table_addr14 + prev_k_End_table_addr14) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr14);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr15 != -1) {
            if ( calAddrtable_addr15( i_Start - prev_i_Start_table_addr15 + prev_i_End_table_addr15, j_Start - prev_j_Start_table_addr15 + prev_j_End_table_addr15, k_Start - prev_k_Start_table_addr15 + prev_k_End_table_addr15) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr15);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr12 = cnt;
                            prev_i_Start_table_addr12 = i_Start;
                            prev_i_End_table_addr12 = i;
                            prev_j_Start_table_addr12 = j_Start;
                            prev_j_End_table_addr12 = j;
                            prev_k_Start_table_addr12 = k_Start;
                            prev_k_End_table_addr12 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr13 = cnt;
                            prev_i_Start_table_addr13 = i_Start;
                            prev_i_End_table_addr13 = i;
                            prev_j_Start_table_addr13 = j_Start;
                            prev_j_End_table_addr13 = j;
                            prev_k_Start_table_addr13 = k_Start;
                            prev_k_End_table_addr13 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr14 = cnt;
                            prev_i_Start_table_addr14 = i_Start;
                            prev_i_End_table_addr14 = i;
                            prev_j_Start_table_addr14 = j_Start;
                            prev_j_End_table_addr14 = j;
                            prev_k_Start_table_addr14 = k_Start;
                            prev_k_End_table_addr14 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr13(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr15 = cnt;
                            prev_i_Start_table_addr15 = i_Start;
                            prev_i_End_table_addr15 = i;
                            prev_j_Start_table_addr15 = j_Start;
                            prev_j_End_table_addr15 = j;
                            prev_k_Start_table_addr15 = k_Start;
                            prev_k_End_table_addr15 = k;
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
void ref_table_addr14() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr12 = -1;
    uint64_t prev_i_Start_table_addr12 = -1;
    uint64_t prev_i_End_table_addr12 = -1;
    uint64_t prev_j_Start_table_addr12 = -1;
    uint64_t prev_j_End_table_addr12 = -1;
    uint64_t prev_k_Start_table_addr12 = -1;
    uint64_t prev_k_End_table_addr12 = -1;
    uint64_t prev_cnt_table_addr13 = -1;
    uint64_t prev_i_Start_table_addr13 = -1;
    uint64_t prev_i_End_table_addr13 = -1;
    uint64_t prev_j_Start_table_addr13 = -1;
    uint64_t prev_j_End_table_addr13 = -1;
    uint64_t prev_k_Start_table_addr13 = -1;
    uint64_t prev_k_End_table_addr13 = -1;
    uint64_t prev_cnt_table_addr14 = -1;
    uint64_t prev_i_Start_table_addr14 = -1;
    uint64_t prev_i_End_table_addr14 = -1;
    uint64_t prev_j_Start_table_addr14 = -1;
    uint64_t prev_j_End_table_addr14 = -1;
    uint64_t prev_k_Start_table_addr14 = -1;
    uint64_t prev_k_End_table_addr14 = -1;
    uint64_t prev_cnt_table_addr15 = -1;
    uint64_t prev_i_Start_table_addr15 = -1;
    uint64_t prev_i_End_table_addr15 = -1;
    uint64_t prev_j_Start_table_addr15 = -1;
    uint64_t prev_j_End_table_addr15 = -1;
    uint64_t prev_k_Start_table_addr15 = -1;
    uint64_t prev_k_End_table_addr15 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        if ( (j_Start - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr12 != -1) {
            if ( calAddrtable_addr12( i_Start - prev_i_Start_table_addr12 + prev_i_End_table_addr12, j_Start - prev_j_Start_table_addr12 + prev_j_End_table_addr12, k_Start - prev_k_Start_table_addr12 + prev_k_End_table_addr12) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr12);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr13 != -1) {
            if ( calAddrtable_addr13( i_Start - prev_i_Start_table_addr13 + prev_i_End_table_addr13, j_Start - prev_j_Start_table_addr13 + prev_j_End_table_addr13, k_Start - prev_k_Start_table_addr13 + prev_k_End_table_addr13) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr13);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr14 != -1) {
            if ( calAddrtable_addr14( i_Start - prev_i_Start_table_addr14 + prev_i_End_table_addr14, j_Start - prev_j_Start_table_addr14 + prev_j_End_table_addr14, k_Start - prev_k_Start_table_addr14 + prev_k_End_table_addr14) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr14);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr15 != -1) {
            if ( calAddrtable_addr15( i_Start - prev_i_Start_table_addr15 + prev_i_End_table_addr15, j_Start - prev_j_Start_table_addr15 + prev_j_End_table_addr15, k_Start - prev_k_Start_table_addr15 + prev_k_End_table_addr15) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr15);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr12 = cnt;
                            prev_i_Start_table_addr12 = i_Start;
                            prev_i_End_table_addr12 = i;
                            prev_j_Start_table_addr12 = j_Start;
                            prev_j_End_table_addr12 = j;
                            prev_k_Start_table_addr12 = k_Start;
                            prev_k_End_table_addr12 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr13 = cnt;
                            prev_i_Start_table_addr13 = i_Start;
                            prev_i_End_table_addr13 = i;
                            prev_j_Start_table_addr13 = j_Start;
                            prev_j_End_table_addr13 = j;
                            prev_k_Start_table_addr13 = k_Start;
                            prev_k_End_table_addr13 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr14 = cnt;
                            prev_i_Start_table_addr14 = i_Start;
                            prev_i_End_table_addr14 = i;
                            prev_j_Start_table_addr14 = j_Start;
                            prev_j_End_table_addr14 = j;
                            prev_k_Start_table_addr14 = k_Start;
                            prev_k_End_table_addr14 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr15 = cnt;
                            prev_i_Start_table_addr15 = i_Start;
                            prev_i_End_table_addr15 = i;
                            prev_j_Start_table_addr15 = j_Start;
                            prev_j_End_table_addr15 = j;
                            prev_k_Start_table_addr15 = k_Start;
                            prev_k_End_table_addr15 = k;
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
void ref_table_addr15() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_table_addr12 = -1;
    uint64_t prev_i_Start_table_addr12 = -1;
    uint64_t prev_i_End_table_addr12 = -1;
    uint64_t prev_j_Start_table_addr12 = -1;
    uint64_t prev_j_End_table_addr12 = -1;
    uint64_t prev_k_Start_table_addr12 = -1;
    uint64_t prev_k_End_table_addr12 = -1;
    uint64_t prev_cnt_table_addr13 = -1;
    uint64_t prev_i_Start_table_addr13 = -1;
    uint64_t prev_i_End_table_addr13 = -1;
    uint64_t prev_j_Start_table_addr13 = -1;
    uint64_t prev_j_End_table_addr13 = -1;
    uint64_t prev_k_Start_table_addr13 = -1;
    uint64_t prev_k_End_table_addr13 = -1;
    uint64_t prev_cnt_table_addr14 = -1;
    uint64_t prev_i_Start_table_addr14 = -1;
    uint64_t prev_i_End_table_addr14 = -1;
    uint64_t prev_j_Start_table_addr14 = -1;
    uint64_t prev_j_End_table_addr14 = -1;
    uint64_t prev_k_Start_table_addr14 = -1;
    uint64_t prev_k_End_table_addr14 = -1;
    uint64_t prev_cnt_table_addr15 = -1;
    uint64_t prev_i_Start_table_addr15 = -1;
    uint64_t prev_i_End_table_addr15 = -1;
    uint64_t prev_j_Start_table_addr15 = -1;
    uint64_t prev_j_End_table_addr15 = -1;
    uint64_t prev_k_Start_table_addr15 = -1;
    uint64_t prev_k_End_table_addr15 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (0 - 1023 + 1) + 1023;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        if ( (j_Start - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_table_addr12 != -1) {
            if ( calAddrtable_addr12( i_Start - prev_i_Start_table_addr12 + prev_i_End_table_addr12, j_Start - prev_j_Start_table_addr12 + prev_j_End_table_addr12, k_Start - prev_k_Start_table_addr12 + prev_k_End_table_addr12) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr12);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr13 != -1) {
            if ( calAddrtable_addr13( i_Start - prev_i_Start_table_addr13 + prev_i_End_table_addr13, j_Start - prev_j_Start_table_addr13 + prev_j_End_table_addr13, k_Start - prev_k_Start_table_addr13 + prev_k_End_table_addr13) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr13);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr14 != -1) {
            if ( calAddrtable_addr14( i_Start - prev_i_Start_table_addr14 + prev_i_End_table_addr14, j_Start - prev_j_Start_table_addr14 + prev_j_End_table_addr14, k_Start - prev_k_Start_table_addr14 + prev_k_End_table_addr14) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr14);
                goto EndSample;
            }
        }
        if ( prev_cnt_table_addr15 != -1) {
            if ( calAddrtable_addr15( i_Start - prev_i_Start_table_addr15 + prev_i_End_table_addr15, j_Start - prev_j_Start_table_addr15 + prev_j_End_table_addr15, k_Start - prev_k_Start_table_addr15 + prev_k_End_table_addr15) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_table_addr15);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i >= 0; i--) {
            {
            int jLB1 = (i + 1);
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr8( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr9( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr12( i, j, k) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr12 = cnt;
                            prev_i_Start_table_addr12 = i_Start;
                            prev_i_End_table_addr12 = i;
                            prev_j_Start_table_addr12 = j_Start;
                            prev_j_End_table_addr12 = j;
                            prev_k_Start_table_addr12 = k_Start;
                            prev_k_End_table_addr12 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr13( i, j, k) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr13 = cnt;
                            prev_i_Start_table_addr13 = i_Start;
                            prev_i_End_table_addr13 = i;
                            prev_j_Start_table_addr13 = j_Start;
                            prev_j_End_table_addr13 = j;
                            prev_k_Start_table_addr13 = k_Start;
                            prev_k_End_table_addr13 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr14 = cnt;
                            prev_i_Start_table_addr14 = i_Start;
                            prev_i_End_table_addr14 = i;
                            prev_j_Start_table_addr14 = j_Start;
                            prev_j_End_table_addr14 = j;
                            prev_k_Start_table_addr14 = k_Start;
                            prev_k_End_table_addr14 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_table_addr15 = cnt;
                            prev_i_Start_table_addr15 = i_Start;
                            prev_i_End_table_addr15 = i;
                            prev_j_Start_table_addr15 = j_Start;
                            prev_j_End_table_addr15 = j;
                            prev_k_Start_table_addr15 = k_Start;
                            prev_k_End_table_addr15 = k;
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
    ref_seq_addr0();
    ref_seq_addr1();
    ref_table_addr0();
    ref_table_addr1();
    ref_table_addr2();
    ref_table_addr3();
    ref_table_addr4();
    ref_table_addr5();
    ref_table_addr6();
    ref_table_addr7();
    ref_table_addr8();
    ref_table_addr9();
    ref_table_addr10();
    ref_table_addr11();
    ref_table_addr12();
    ref_table_addr13();
    ref_table_addr14();
    ref_table_addr15();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
nussinov */ 
