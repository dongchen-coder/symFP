
 /* Start to analysis array index
Array index info: Total number of references: 18
table.addr (((i + 1) * 1024) + j)
table.addr ((i * 1024) + j)
table.addr (((i * 1024) + j) - 1)
table.addr ((i * 1024) + j)
table.addr ((((i + 1) * 1024) + j) - 1)
seq.addr i
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((((i + 1) * 1024) + j) - 1)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
seq.addr j
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + k)
table.addr (((k + 1) * 1024) + j)

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
/* table_addr ((i * 1024) + j) 0 */
int calAddrtable_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr (((i * 1024) + j) - 1) 1 */
int calAddrtable_addr1( int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + j) 2 */
int calAddrtable_addr2( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + j) 3 */
int calAddrtable_addr3( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr (((i + 1) * 1024) + j) 4 */
int calAddrtable_addr4( int i, int j) {
    int result = ((((i + 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + j) 5 */
int calAddrtable_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + j) 6 */
int calAddrtable_addr6( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr ((((i + 1) * 1024) + j) - 1) 7 */
int calAddrtable_addr7( int i, int j) {
    int result = (((((i + 1) * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* seq_addr i 8 */
int calAddrseq_addr8( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* seq_addr j 9 */
int calAddrseq_addr9( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + j) 10 */
int calAddrtable_addr10( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + j) 11 */
int calAddrtable_addr11( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr ((((i + 1) * 1024) + j) - 1) 12 */
int calAddrtable_addr12( int i, int j) {
    int result = (((((i + 1) * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + j) 13 */
int calAddrtable_addr13( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + j) 14 */
int calAddrtable_addr14( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + k) 15 */
int calAddrtable_addr15( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* table_addr (((k + 1) * 1024) + j) 16 */
int calAddrtable_addr16( int i, int j, int k) {
    int result = ((((k + 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* table_addr ((i * 1024) + j) 17 */
int calAddrtable_addr17( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_table_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr4(i_Start, j_Start)) {
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
void ref_table_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr0(i_Start, j_Start)) {
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
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr1(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr1(i_Start, j_Start)) {
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
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr6(i_Start, j_Start)) {
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
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr7(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr7(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr7(i_Start, j_Start)) {
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
void ref_seq_addr8() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                    if ( calAddrseq_addr8( i, j) == calAddrseq_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrseq_addr9( i, j) == calAddrseq_addr8(i_Start, j_Start)) {
                        rtHistoCal(cnt);
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
void ref_table_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr2(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr2(i_Start, j_Start)) {
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
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr3(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr3(i_Start, j_Start)) {
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
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr10(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr10(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr10(i_Start, j_Start)) {
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
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr11(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr11(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr11(i_Start, j_Start)) {
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
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr12(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr12(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr12(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr12(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr12(i_Start, j_Start)) {
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
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr5(i_Start, j_Start)) {
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
void ref_table_addr17() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        if ( (j_Start - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr17(i_Start, j_Start, k_Start)) {
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
void ref_seq_addr9() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                    if ( calAddrseq_addr8( i, j) == calAddrseq_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrseq_addr9( i, j) == calAddrseq_addr9(i_Start, j_Start)) {
                        rtHistoCal(cnt);
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
void ref_table_addr13() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr13(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB2 = (i + 1);
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr13(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr13(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr13(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr13(i_Start, j_Start)) {
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
void ref_table_addr14() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        if ( (j_Start - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr14(i_Start, j_Start, k_Start)) {
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
void ref_table_addr15() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        if ( (j_Start - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr15(i_Start, j_Start, k_Start)) {
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
void ref_table_addr16() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        if ( (1024 - (i_Start + 1)) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - (i_Start + 1)) + (i_Start + 1);
        if ( (j_Start - (i_Start + 1)) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - (i_Start + 1)) + (i_Start + 1);
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

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
                    if ( calAddrtable_addr0( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr1( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr2( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr3( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr4( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr5( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr6( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr7( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr10( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr11( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr12( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtable_addr13( i, j) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
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
                        if ( calAddrtable_addr14( i, j, k) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr15( i, j, k) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr16( i, j, k) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtable_addr17( i, j, k) == calAddrtable_addr16(i_Start, j_Start, k_Start)) {
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
int main() {
    ref_table_addr4();
    ref_table_addr0();
    ref_table_addr1();
    ref_table_addr6();
    ref_table_addr7();
    ref_seq_addr8();
    ref_table_addr2();
    ref_table_addr3();
    ref_table_addr10();
    ref_table_addr11();
    ref_table_addr12();
    ref_table_addr5();
    ref_table_addr17();
    ref_seq_addr9();
    ref_table_addr13();
    ref_table_addr14();
    ref_table_addr15();
    ref_table_addr16();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: nussinov */ 
