
 /* Start to analysis array index
Array index info: Total number of references: 22
A.addr ((((i * 256) * 256) + ((j + 1) * 256)) + k)
A.addr (((((i + 1) * 256) * 256) + (j * 256)) + k)
A.addr ((((i * 256) * 256) + ((j - 1) * 256)) + k)
A.addr ((((i * 256) * 256) + (j * 256)) + k)
A.addr (((((i - 1) * 256) * 256) + (j * 256)) + k)
A.addr ((((i * 256) * 256) + (j * 256)) + k)
A.addr (((((i * 256) * 256) + (j * 256)) + k) - 1)
A.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr (((((i + 1) * 256) * 256) + (j * 256)) + k)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr (((((i - 1) * 256) * 256) + (j * 256)) + k)
A.addr (((((i * 256) * 256) + (j * 256)) + k) + 1)
B.addr ((((i * 256) * 256) + ((j + 1) * 256)) + k)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr ((((i * 256) * 256) + ((j - 1) * 256)) + k)
B.addr (((((i * 256) * 256) + (j * 256)) + k) + 1)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr (((((i * 256) * 256) + (j * 256)) + k) - 1)
A.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
A.addr ((((i * 256) * 256) + (j * 256)) + k)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %B
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (1, 10)
--Loop inc: (t + 1)
--Loop predicate: <=
----i
----Loop Bound: (1, 255)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (1, 255)
------Loop inc: (j + 1)
------Loop predicate: <
--------k
--------Loop Bound: (1, 255)
--------Loop inc: (k + 1)
--------Loop predicate: <
----------array access A.addr (((((i + 1) * 256) * 256) + (j * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access A.addr (((((i - 1) * 256) * 256) + (j * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + ((j + 1) * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + ((j - 1) * 256)) + k)
----------array access A.addr (((((i * 256) * 256) + (j * 256)) + k) + 1)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access A.addr (((((i * 256) * 256) + (j * 256)) + k) - 1)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----i
----Loop Bound: (1, 255)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (1, 255)
------Loop inc: (j + 1)
------Loop predicate: <
--------k
--------Loop Bound: (1, 255)
--------Loop inc: (k + 1)
--------Loop predicate: <
----------array access B.addr (((((i + 1) * 256) * 256) + (j * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access B.addr (((((i - 1) * 256) * 256) + (j * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + ((j + 1) * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + ((j - 1) * 256)) + k)
----------array access B.addr (((((i * 256) * 256) + (j * 256)) + k) + 1)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access B.addr (((((i * 256) * 256) + (j * 256)) + k) - 1)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 0
------Sample number: 0
--------Sample number: 0
----------Sample number: 1
------Sample number: 0
--------Sample number: 0
----------Sample number: 1
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
/* A_addr (((((i + 1) * 256) * 256) + (j * 256)) + k) 0 */
int calAddrA_addr0( int t, int i, int j, int k) {
    int result = ((((((i + 1) * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 1 */
int calAddrA_addr1( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr (((((i - 1) * 256) * 256) + (j * 256)) + k) 2 */
int calAddrA_addr2( int t, int i, int j, int k) {
    int result = ((((((i - 1) * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + ((j + 1) * 256)) + k) 3 */
int calAddrA_addr3( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + ((j + 1) * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 4 */
int calAddrA_addr4( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + ((j - 1) * 256)) + k) 5 */
int calAddrA_addr5( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + ((j - 1) * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr (((((i * 256) * 256) + (j * 256)) + k) + 1) 6 */
int calAddrA_addr6( int t, int i, int j, int k) {
    int result = ((((((i * 256) * 256) + (j * 256)) + k) + 1)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 7 */
int calAddrA_addr7( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr (((((i * 256) * 256) + (j * 256)) + k) - 1) 8 */
int calAddrA_addr8( int t, int i, int j, int k) {
    int result = ((((((i * 256) * 256) + (j * 256)) + k) - 1)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 9 */
int calAddrA_addr9( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 10 */
int calAddrB_addr10( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr (((((i + 1) * 256) * 256) + (j * 256)) + k) 11 */
int calAddrB_addr11( int t, int i, int j, int k) {
    int result = ((((((i + 1) * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 12 */
int calAddrB_addr12( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr (((((i - 1) * 256) * 256) + (j * 256)) + k) 13 */
int calAddrB_addr13( int t, int i, int j, int k) {
    int result = ((((((i - 1) * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + ((j + 1) * 256)) + k) 14 */
int calAddrB_addr14( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + ((j + 1) * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 15 */
int calAddrB_addr15( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + ((j - 1) * 256)) + k) 16 */
int calAddrB_addr16( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + ((j - 1) * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr (((((i * 256) * 256) + (j * 256)) + k) + 1) 17 */
int calAddrB_addr17( int t, int i, int j, int k) {
    int result = ((((((i * 256) * 256) + (j * 256)) + k) + 1)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 18 */
int calAddrB_addr18( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr (((((i * 256) * 256) + (j * 256)) + k) - 1) 19 */
int calAddrB_addr19( int t, int i, int j, int k) {
    int result = ((((((i * 256) * 256) + (j * 256)) + k) - 1)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 20 */
int calAddrB_addr20( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 21 */
int calAddrA_addr21( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
void ref_A_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
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
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
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
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr7() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr8() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr11() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr11(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr12() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr12(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr13() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr13(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr14() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr14(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr15() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr15(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr16() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr16(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr17() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr17(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr18() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr18(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr19() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        cntStart = true;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr19(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr9() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr10() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_B_addr20() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrB_addr10( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr11( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr12( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr13( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr14( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr15( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr16( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr17( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr18( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr19( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr20( t, i, j, k) == calAddrB_addr20(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr21() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 255; i++) {
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 255; j++) {
                    {
                    int kLB3 = 1;
                    for ( int k = kLB3; k < 255; k++) {
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
                                rtHistoCal(cnt);
                                goto EndSample;
                            }
                        }
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
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
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if (cntStart == true) cnt++;
                        if (cntStart == true) cnt++;
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
                            if ( calAddrA_addr21( t, i, j, k) == calAddrA_addr21(t_Start, i_Start, j_Start, k_Start)) {
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
        }
        }
EndSample:
        s++;
        }
}
int main() {
    ref_A_addr3();
    ref_A_addr0();
    ref_A_addr5();
    ref_A_addr1();
    ref_A_addr2();
    ref_A_addr7();
    ref_A_addr8();
    ref_A_addr4();
    ref_B_addr11();
    ref_B_addr12();
    ref_B_addr13();
    ref_A_addr6();
    ref_B_addr14();
    ref_B_addr15();
    ref_B_addr16();
    ref_B_addr17();
    ref_B_addr18();
    ref_B_addr19();
    ref_A_addr9();
    ref_B_addr10();
    ref_B_addr20();
    ref_A_addr21();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: heat_3d */ 
