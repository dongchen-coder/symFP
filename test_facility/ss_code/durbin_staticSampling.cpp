
 /* Start to analysis array index
Array index info
r.addr ((k - i) - 1)
y.addr i
r.addr k
y.addr i
y.addr ((k - i) - 1)
z.addr i
z.addr i
y.addr i
y.addr k

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %y
double* %r
double* %z

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--k
--Loop Bound: (1, 1024)
--Loop inc: (k + 1)
--Loop predicate: <
----i
----Loop Bound: (0, k)
----Loop inc: (i + 1)
----Loop predicate: <
------array access r.addr ((k - i) - 1)
------array access y.addr i
----array access r.addr k
----i
----Loop Bound: (0, k)
----Loop inc: (i + 1)
----Loop predicate: <
------array access y.addr i
------array access y.addr ((k - i) - 1)
------array access z.addr i
----i
----Loop Bound: (0, k)
----Loop inc: (i + 1)
----Loop predicate: <
------array access z.addr i
------array access y.addr i
----array access y.addr k

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 1 0 
Dump stride: 1 1 
init counter: 1 0 
Dump stride: 1 1 
init counter: 1 0 
Dump stride: 1 1 
Dump tree:
----Sample number: 10
------Sample number: 52
------Sample number: 52
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
int calAddrr_addr0( int k, int i) {
    int result = (((k - i) - 1)) * 8 / 64;
    return result;
}
int calAddry_addr0( int k, int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrr_addr1( int k) {
    int result = (k) * 8 / 64;
    return result;
}
int calAddry_addr1( int k, int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddry_addr2( int k, int i) {
    int result = (((k - i) - 1)) * 8 / 64;
    return result;
}
int calAddrz_addr0( int k, int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddrz_addr1( int k, int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddry_addr3( int k, int i) {
    int result = (i) * 8 / 64;
    return result;
}
int calAddry_addr4( int k) {
    int result = (k) * 8 / 64;
    return result;
}
void ref_r_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int k_Start = rand() % (1024 - 1) + 1;
        if ( (k_Start - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (k_Start - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < k; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrr_addr0( k, i) == calAddrr_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrr_addr1( k) == calAddrr_addr0(k_Start, i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
        }
        }
EndSample:
        s++;
        }
}
void ref_r_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int k_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < k; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrr_addr0( k, i) == calAddrr_addr1(k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrr_addr1( k) == calAddrr_addr1(k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
        }
        }
EndSample:
        s++;
        }
}
void ref_y_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int k_Start = rand() % (1024 - 1) + 1;
        if ( (k_Start - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (k_Start - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr0( k, i) == calAddry_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < k; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr1( k, i) == calAddry_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr2( k, i) == calAddry_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr3( k, i) == calAddry_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddry_addr4( k) == calAddry_addr0(k_Start, i_Start)) {
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
void ref_y_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int k_Start = rand() % (1024 - 1) + 1;
        if ( (k_Start - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (k_Start - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr0( k, i) == calAddry_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            if ( k == k_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < k; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr1( k, i) == calAddry_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr2( k, i) == calAddry_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr3( k, i) == calAddry_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddry_addr4( k) == calAddry_addr1(k_Start, i_Start)) {
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
void ref_y_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int k_Start = rand() % (1024 - 1) + 1;
        if ( (k_Start - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (k_Start - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr0( k, i) == calAddry_addr2(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            if ( k == k_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < k; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr1( k, i) == calAddry_addr2(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr2( k, i) == calAddry_addr2(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr3( k, i) == calAddry_addr2(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddry_addr4( k) == calAddry_addr2(k_Start, i_Start)) {
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
void ref_y_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int k_Start = rand() % (1024 - 1) + 1;
        if ( (k_Start - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (k_Start - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr0( k, i) == calAddry_addr3(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < k; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr1( k, i) == calAddry_addr3(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr2( k, i) == calAddry_addr3(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB3 = 0;
            if ( k == k_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr3( k, i) == calAddry_addr3(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddry_addr4( k) == calAddry_addr3(k_Start, i_Start)) {
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
void ref_y_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int k_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr0( k, i) == calAddry_addr4(k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < k; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr1( k, i) == calAddry_addr4(k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr2( k, i) == calAddry_addr4(k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr3( k, i) == calAddry_addr4(k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddry_addr4( k) == calAddry_addr4(k_Start)) {
                    rtHistoCal(cnt);
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
void ref_z_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int k_Start = rand() % (1024 - 1) + 1;
        if ( (k_Start - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (k_Start - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            if ( k == k_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrz_addr0( k, i) == calAddrz_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            {
            int iLB3 = 0;
            for ( int i = iLB3; i < k; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrz_addr1( k, i) == calAddrz_addr0(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
        }
        }
EndSample:
        s++;
        }
}
void ref_z_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int k_Start = rand() % (1024 - 1) + 1;
        if ( (k_Start - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (k_Start - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            for ( int i = iLB1; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int iLB2 = 0;
            for ( int i = iLB2; i < k; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrz_addr0( k, i) == calAddrz_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int iLB3 = 0;
            if ( k == k_Start ) {
                iLB3 = i_Start;
            }
            for ( int i = iLB3; i < k; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrz_addr1( k, i) == calAddrz_addr1(k_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
        }
        }
EndSample:
        s++;
        }
}
int main() {
    ref_r_addr0();
    ref_r_addr1();
    ref_y_addr0();
    ref_y_addr1();
    ref_y_addr2();
    ref_y_addr3();
    ref_y_addr4();
    ref_z_addr0();
    ref_z_addr1();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: durbin */ 
