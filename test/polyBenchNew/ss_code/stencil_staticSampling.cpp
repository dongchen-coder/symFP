
 /* Start to analysis array index
Array index info
b.addr ((1026 * i) + j)
a.addr ((1026 * i) + j)
a.addr (((1026 * i) + j) + 1)
a.addr (((1026 * i) + j) - 1)
a.addr ((1026 * (i - 1)) + j)
a.addr ((1026 * (i + 1)) + j)

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
----j
----Loop Bound: (1, 1025)
------array access a.addr ((1026 * i) + j)
------array access a.addr (((1026 * i) + j) + 1)
------array access a.addr (((1026 * i) + j) - 1)
------array access a.addr ((1026 * (i - 1)) + j)
------array access a.addr ((1026 * (i + 1)) + j)
------array access b.addr ((1026 * i) + j)

Finish analysis loops */ 
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
    for (uint64_t c = 0; c <= max_RT; c++) {
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
    int result = (((1026 * i) + j)) * 8 / 32;
    return result;
}
int calAddra_addr1( int i, int j) {
    int result = ((((1026 * i) + j) + 1)) * 8 / 32;
    return result;
}
int calAddra_addr2( int i, int j) {
    int result = ((((1026 * i) + j) - 1)) * 8 / 32;
    return result;
}
int calAddra_addr3( int i, int j) {
    int result = (((1026 * (i - 1)) + j)) * 8 / 32;
    return result;
}
int calAddra_addr4( int i, int j) {
    int result = (((1026 * (i + 1)) + j)) * 8 / 32;
    return result;
}
int calAddrb_addr0( int i, int j) {
    int result = (((1026 * i) + j)) * 8 / 32;
    return result;
}
void updateCoefficient(float* a, int** x, int length, int* y) {
    double** x_tmp = new double*[length];
    for (int i = 0; i < length; i++) {
        x_tmp[i] = new double[length];
        for (int j = 0; j < length; j++) {
            x_tmp[i][j] = (double) x[i][j];
        }
    }
    double* y_tmp = new double[length];
    for (int i = 0; i < length; i++) {
        y_tmp[i] = y[i];
    }
    for (int i = 0; i < length; i++) {
        double maxEl = abs(x_tmp[i][i]);
        int maxRow = i;
        for (int k=i+1; k< length; k++) {
            if (abs(x_tmp[k][i]) > maxEl) {
                maxEl = abs(x_tmp[k][i]);
                maxRow = k;
            }
        }
        for (int k=i; k< length;k++) {
            double tmp = x_tmp[maxRow][k];
            x_tmp[maxRow][k] = x_tmp[i][k];
            x_tmp[i][k] = tmp;
        }
        double tmp = y_tmp[maxRow];
        y_tmp[maxRow] = y_tmp[i];
        y_tmp[i] = tmp;
        for (int k=i+1; k< length; k++) {
            double c = -x_tmp[i][k]/x_tmp[i][i];
            for (int j=i; j< length; j++) {
                if (i==j) {
                    x_tmp[j][k] = 0;
                } else {
                    x_tmp[j][k] += c * x_tmp[j][i];
                }
            }
            y_tmp[k] += c * y_tmp[i];
        }
    }
    for (int i=length-1; i>=0; i--) {
        a[i] = y_tmp[i]/x_tmp[i][i];
        for (int k=i-1; k>=0; k--) {
            y_tmp[k] -= x_tmp[i][k] * a[i];
        }
    }
}
int calWithCoefficient(float* a, int *x, int length) {
    float tmp = 0;
    for (int i = 0; i < length; i++) {
        tmp += a[i] * x[i];
    }
    tmp += a[length];
    return (int) tmp;
}
void ref_a_addr0() {
    /* Generating profile counter */
    uint64_t pred_num = 0;
    uint64_t pred_sl_num = 0;
    uint64_t pred_dl_num = 0;
    uint64_t sample_num = 0;
    uint64_t reuse_in_same_loop = 0;
    uint64_t reuse_in_diff_loops = 0;
    uint64_t no_reuse_num = 0;
    uint64_t mis_pred_num = 0;
    uint64_t pred = 0;
    uint64_t actural = 0;
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_a_addr0 = -1;
    uint64_t prev_i_Start_a_addr0 = -1;
    uint64_t prev_i_End_a_addr0 = -1;
    uint64_t prev_j_Start_a_addr0 = -1;
    uint64_t prev_j_End_a_addr0 = -1;
    uint64_t prev_cnt_a_addr1 = -1;
    uint64_t prev_i_Start_a_addr1 = -1;
    uint64_t prev_i_End_a_addr1 = -1;
    uint64_t prev_j_Start_a_addr1 = -1;
    uint64_t prev_j_End_a_addr1 = -1;
    uint64_t prev_cnt_a_addr2 = -1;
    uint64_t prev_i_Start_a_addr2 = -1;
    uint64_t prev_i_End_a_addr2 = -1;
    uint64_t prev_j_Start_a_addr2 = -1;
    uint64_t prev_j_End_a_addr2 = -1;
    uint64_t prev_cnt_a_addr3 = -1;
    uint64_t prev_i_Start_a_addr3 = -1;
    uint64_t prev_i_End_a_addr3 = -1;
    uint64_t prev_j_Start_a_addr3 = -1;
    uint64_t prev_j_End_a_addr3 = -1;
    uint64_t prev_cnt_a_addr4 = -1;
    uint64_t prev_i_Start_a_addr4 = -1;
    uint64_t prev_i_End_a_addr4 = -1;
    uint64_t prev_j_Start_a_addr4 = -1;
    uint64_t prev_j_End_a_addr4 = -1;
    /* Generating search reuse init code (different loops) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1025 - 1) + 1;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1025 - 1) + 1;
            j_Start = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_a_addr0 != -1) {
            if ( calAddra_addr0( i_Start - prev_i_Start_a_addr0 + prev_i_End_a_addr0, j_Start - prev_j_Start_a_addr0 + prev_j_End_a_addr0) == calAddra_addr0(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr0;
                rtHistoCal(prev_cnt_a_addr0);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr1 != -1) {
            if ( calAddra_addr1( i_Start - prev_i_Start_a_addr1 + prev_i_End_a_addr1, j_Start - prev_j_Start_a_addr1 + prev_j_End_a_addr1) == calAddra_addr0(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr1;
                rtHistoCal(prev_cnt_a_addr1);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr2 != -1) {
            if ( calAddra_addr2( i_Start - prev_i_Start_a_addr2 + prev_i_End_a_addr2, j_Start - prev_j_Start_a_addr2 + prev_j_End_a_addr2) == calAddra_addr0(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr2;
                rtHistoCal(prev_cnt_a_addr2);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr3 != -1) {
            if ( calAddra_addr3( i_Start - prev_i_Start_a_addr3 + prev_i_End_a_addr3, j_Start - prev_j_Start_a_addr3 + prev_j_End_a_addr3) == calAddra_addr0(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr3;
                rtHistoCal(prev_cnt_a_addr3);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr4 != -1) {
            if ( calAddra_addr4( i_Start - prev_i_Start_a_addr4 + prev_i_End_a_addr4, j_Start - prev_j_Start_a_addr4 + prev_j_End_a_addr4) == calAddra_addr0(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr4;
                rtHistoCal(prev_cnt_a_addr4);
                goto PROFILE;
            }
        }
        /* Generating search reuse body code (use reuse are in different loop) */
        /* Finished search reuse */
PROFILE:
sample_num++;
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
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr0 = cnt;
                        prev_i_Start_a_addr0 = i_Start;
                        prev_i_End_a_addr0 = i;
                        prev_j_Start_a_addr0 = j_Start;
                        prev_j_End_a_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr1 = cnt;
                        prev_i_Start_a_addr1 = i_Start;
                        prev_i_End_a_addr1 = i;
                        prev_j_Start_a_addr1 = j_Start;
                        prev_j_End_a_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr2 = cnt;
                        prev_i_Start_a_addr2 = i_Start;
                        prev_i_End_a_addr2 = i;
                        prev_j_Start_a_addr2 = j_Start;
                        prev_j_End_a_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr3 = cnt;
                        prev_i_Start_a_addr3 = i_Start;
                        prev_i_End_a_addr3 = i;
                        prev_j_Start_a_addr3 = j_Start;
                        prev_j_End_a_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr4 = cnt;
                        prev_i_Start_a_addr4 = i_Start;
                        prev_i_End_a_addr4 = i;
                        prev_j_Start_a_addr4 = j_Start;
                        prev_j_End_a_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        no_reuse_num++;
        if (actural != 0 && actural != pred) mis_pred_num++;
        actural = 0;
EndSample:
        s++;
        }
std::cout << "SN " << sample_num << " ";
std::cout << "NSR " << sample_num - no_reuse_num << " ";
std::cout << "RSL " << reuse_in_same_loop << " ";
std::cout << "RDL " << reuse_in_diff_loops << " ";
std::cout << "PN " << pred_num << " ";
std::cout << "PSL " << pred_sl_num << " ";
std::cout << "PDL " << pred_dl_num << " ";
std::cout << "MPN " << mis_pred_num << std::endl;
}
void ref_a_addr1() {
    /* Generating profile counter */
    uint64_t pred_num = 0;
    uint64_t pred_sl_num = 0;
    uint64_t pred_dl_num = 0;
    uint64_t sample_num = 0;
    uint64_t reuse_in_same_loop = 0;
    uint64_t reuse_in_diff_loops = 0;
    uint64_t no_reuse_num = 0;
    uint64_t mis_pred_num = 0;
    uint64_t pred = 0;
    uint64_t actural = 0;
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_a_addr0 = -1;
    uint64_t prev_i_Start_a_addr0 = -1;
    uint64_t prev_i_End_a_addr0 = -1;
    uint64_t prev_j_Start_a_addr0 = -1;
    uint64_t prev_j_End_a_addr0 = -1;
    uint64_t prev_cnt_a_addr1 = -1;
    uint64_t prev_i_Start_a_addr1 = -1;
    uint64_t prev_i_End_a_addr1 = -1;
    uint64_t prev_j_Start_a_addr1 = -1;
    uint64_t prev_j_End_a_addr1 = -1;
    uint64_t prev_cnt_a_addr2 = -1;
    uint64_t prev_i_Start_a_addr2 = -1;
    uint64_t prev_i_End_a_addr2 = -1;
    uint64_t prev_j_Start_a_addr2 = -1;
    uint64_t prev_j_End_a_addr2 = -1;
    uint64_t prev_cnt_a_addr3 = -1;
    uint64_t prev_i_Start_a_addr3 = -1;
    uint64_t prev_i_End_a_addr3 = -1;
    uint64_t prev_j_Start_a_addr3 = -1;
    uint64_t prev_j_End_a_addr3 = -1;
    uint64_t prev_cnt_a_addr4 = -1;
    uint64_t prev_i_Start_a_addr4 = -1;
    uint64_t prev_i_End_a_addr4 = -1;
    uint64_t prev_j_Start_a_addr4 = -1;
    uint64_t prev_j_End_a_addr4 = -1;
    /* Generating search reuse init code (different loops) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1025 - 1) + 1;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1025 - 1) + 1;
            j_Start = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_a_addr0 != -1) {
            if ( calAddra_addr0( i_Start - prev_i_Start_a_addr0 + prev_i_End_a_addr0, j_Start - prev_j_Start_a_addr0 + prev_j_End_a_addr0) == calAddra_addr1(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr0;
                rtHistoCal(prev_cnt_a_addr0);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr1 != -1) {
            if ( calAddra_addr1( i_Start - prev_i_Start_a_addr1 + prev_i_End_a_addr1, j_Start - prev_j_Start_a_addr1 + prev_j_End_a_addr1) == calAddra_addr1(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr1;
                rtHistoCal(prev_cnt_a_addr1);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr2 != -1) {
            if ( calAddra_addr2( i_Start - prev_i_Start_a_addr2 + prev_i_End_a_addr2, j_Start - prev_j_Start_a_addr2 + prev_j_End_a_addr2) == calAddra_addr1(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr2;
                rtHistoCal(prev_cnt_a_addr2);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr3 != -1) {
            if ( calAddra_addr3( i_Start - prev_i_Start_a_addr3 + prev_i_End_a_addr3, j_Start - prev_j_Start_a_addr3 + prev_j_End_a_addr3) == calAddra_addr1(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr3;
                rtHistoCal(prev_cnt_a_addr3);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr4 != -1) {
            if ( calAddra_addr4( i_Start - prev_i_Start_a_addr4 + prev_i_End_a_addr4, j_Start - prev_j_Start_a_addr4 + prev_j_End_a_addr4) == calAddra_addr1(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr4;
                rtHistoCal(prev_cnt_a_addr4);
                goto PROFILE;
            }
        }
        /* Generating search reuse body code (use reuse are in different loop) */
        /* Finished search reuse */
PROFILE:
sample_num++;
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
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr0 = cnt;
                        prev_i_Start_a_addr0 = i_Start;
                        prev_i_End_a_addr0 = i;
                        prev_j_Start_a_addr0 = j_Start;
                        prev_j_End_a_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr1 = cnt;
                        prev_i_Start_a_addr1 = i_Start;
                        prev_i_End_a_addr1 = i;
                        prev_j_Start_a_addr1 = j_Start;
                        prev_j_End_a_addr1 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr2 = cnt;
                        prev_i_Start_a_addr2 = i_Start;
                        prev_i_End_a_addr2 = i;
                        prev_j_Start_a_addr2 = j_Start;
                        prev_j_End_a_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr3 = cnt;
                        prev_i_Start_a_addr3 = i_Start;
                        prev_i_End_a_addr3 = i;
                        prev_j_Start_a_addr3 = j_Start;
                        prev_j_End_a_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr4 = cnt;
                        prev_i_Start_a_addr4 = i_Start;
                        prev_i_End_a_addr4 = i;
                        prev_j_Start_a_addr4 = j_Start;
                        prev_j_End_a_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        no_reuse_num++;
        if (actural != 0 && actural != pred) mis_pred_num++;
        actural = 0;
EndSample:
        s++;
        }
std::cout << "SN " << sample_num << " ";
std::cout << "NSR " << sample_num - no_reuse_num << " ";
std::cout << "RSL " << reuse_in_same_loop << " ";
std::cout << "RDL " << reuse_in_diff_loops << " ";
std::cout << "PN " << pred_num << " ";
std::cout << "PSL " << pred_sl_num << " ";
std::cout << "PDL " << pred_dl_num << " ";
std::cout << "MPN " << mis_pred_num << std::endl;
}
void ref_a_addr2() {
    /* Generating profile counter */
    uint64_t pred_num = 0;
    uint64_t pred_sl_num = 0;
    uint64_t pred_dl_num = 0;
    uint64_t sample_num = 0;
    uint64_t reuse_in_same_loop = 0;
    uint64_t reuse_in_diff_loops = 0;
    uint64_t no_reuse_num = 0;
    uint64_t mis_pred_num = 0;
    uint64_t pred = 0;
    uint64_t actural = 0;
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_a_addr0 = -1;
    uint64_t prev_i_Start_a_addr0 = -1;
    uint64_t prev_i_End_a_addr0 = -1;
    uint64_t prev_j_Start_a_addr0 = -1;
    uint64_t prev_j_End_a_addr0 = -1;
    uint64_t prev_cnt_a_addr1 = -1;
    uint64_t prev_i_Start_a_addr1 = -1;
    uint64_t prev_i_End_a_addr1 = -1;
    uint64_t prev_j_Start_a_addr1 = -1;
    uint64_t prev_j_End_a_addr1 = -1;
    uint64_t prev_cnt_a_addr2 = -1;
    uint64_t prev_i_Start_a_addr2 = -1;
    uint64_t prev_i_End_a_addr2 = -1;
    uint64_t prev_j_Start_a_addr2 = -1;
    uint64_t prev_j_End_a_addr2 = -1;
    uint64_t prev_cnt_a_addr3 = -1;
    uint64_t prev_i_Start_a_addr3 = -1;
    uint64_t prev_i_End_a_addr3 = -1;
    uint64_t prev_j_Start_a_addr3 = -1;
    uint64_t prev_j_End_a_addr3 = -1;
    uint64_t prev_cnt_a_addr4 = -1;
    uint64_t prev_i_Start_a_addr4 = -1;
    uint64_t prev_i_End_a_addr4 = -1;
    uint64_t prev_j_Start_a_addr4 = -1;
    uint64_t prev_j_End_a_addr4 = -1;
    /* Generating search reuse init code (different loops) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1025 - 1) + 1;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1025 - 1) + 1;
            j_Start = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_a_addr0 != -1) {
            if ( calAddra_addr0( i_Start - prev_i_Start_a_addr0 + prev_i_End_a_addr0, j_Start - prev_j_Start_a_addr0 + prev_j_End_a_addr0) == calAddra_addr2(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr0;
                rtHistoCal(prev_cnt_a_addr0);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr1 != -1) {
            if ( calAddra_addr1( i_Start - prev_i_Start_a_addr1 + prev_i_End_a_addr1, j_Start - prev_j_Start_a_addr1 + prev_j_End_a_addr1) == calAddra_addr2(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr1;
                rtHistoCal(prev_cnt_a_addr1);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr2 != -1) {
            if ( calAddra_addr2( i_Start - prev_i_Start_a_addr2 + prev_i_End_a_addr2, j_Start - prev_j_Start_a_addr2 + prev_j_End_a_addr2) == calAddra_addr2(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr2;
                rtHistoCal(prev_cnt_a_addr2);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr3 != -1) {
            if ( calAddra_addr3( i_Start - prev_i_Start_a_addr3 + prev_i_End_a_addr3, j_Start - prev_j_Start_a_addr3 + prev_j_End_a_addr3) == calAddra_addr2(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr3;
                rtHistoCal(prev_cnt_a_addr3);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr4 != -1) {
            if ( calAddra_addr4( i_Start - prev_i_Start_a_addr4 + prev_i_End_a_addr4, j_Start - prev_j_Start_a_addr4 + prev_j_End_a_addr4) == calAddra_addr2(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr4;
                rtHistoCal(prev_cnt_a_addr4);
                goto PROFILE;
            }
        }
        /* Generating search reuse body code (use reuse are in different loop) */
        /* Finished search reuse */
PROFILE:
sample_num++;
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
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr0 = cnt;
                        prev_i_Start_a_addr0 = i_Start;
                        prev_i_End_a_addr0 = i;
                        prev_j_Start_a_addr0 = j_Start;
                        prev_j_End_a_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr1 = cnt;
                        prev_i_Start_a_addr1 = i_Start;
                        prev_i_End_a_addr1 = i;
                        prev_j_Start_a_addr1 = j_Start;
                        prev_j_End_a_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr2 = cnt;
                        prev_i_Start_a_addr2 = i_Start;
                        prev_i_End_a_addr2 = i;
                        prev_j_Start_a_addr2 = j_Start;
                        prev_j_End_a_addr2 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr3 = cnt;
                        prev_i_Start_a_addr3 = i_Start;
                        prev_i_End_a_addr3 = i;
                        prev_j_Start_a_addr3 = j_Start;
                        prev_j_End_a_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr4 = cnt;
                        prev_i_Start_a_addr4 = i_Start;
                        prev_i_End_a_addr4 = i;
                        prev_j_Start_a_addr4 = j_Start;
                        prev_j_End_a_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        no_reuse_num++;
        if (actural != 0 && actural != pred) mis_pred_num++;
        actural = 0;
EndSample:
        s++;
        }
std::cout << "SN " << sample_num << " ";
std::cout << "NSR " << sample_num - no_reuse_num << " ";
std::cout << "RSL " << reuse_in_same_loop << " ";
std::cout << "RDL " << reuse_in_diff_loops << " ";
std::cout << "PN " << pred_num << " ";
std::cout << "PSL " << pred_sl_num << " ";
std::cout << "PDL " << pred_dl_num << " ";
std::cout << "MPN " << mis_pred_num << std::endl;
}
void ref_a_addr3() {
    /* Generating profile counter */
    uint64_t pred_num = 0;
    uint64_t pred_sl_num = 0;
    uint64_t pred_dl_num = 0;
    uint64_t sample_num = 0;
    uint64_t reuse_in_same_loop = 0;
    uint64_t reuse_in_diff_loops = 0;
    uint64_t no_reuse_num = 0;
    uint64_t mis_pred_num = 0;
    uint64_t pred = 0;
    uint64_t actural = 0;
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_a_addr0 = -1;
    uint64_t prev_i_Start_a_addr0 = -1;
    uint64_t prev_i_End_a_addr0 = -1;
    uint64_t prev_j_Start_a_addr0 = -1;
    uint64_t prev_j_End_a_addr0 = -1;
    uint64_t prev_cnt_a_addr1 = -1;
    uint64_t prev_i_Start_a_addr1 = -1;
    uint64_t prev_i_End_a_addr1 = -1;
    uint64_t prev_j_Start_a_addr1 = -1;
    uint64_t prev_j_End_a_addr1 = -1;
    uint64_t prev_cnt_a_addr2 = -1;
    uint64_t prev_i_Start_a_addr2 = -1;
    uint64_t prev_i_End_a_addr2 = -1;
    uint64_t prev_j_Start_a_addr2 = -1;
    uint64_t prev_j_End_a_addr2 = -1;
    uint64_t prev_cnt_a_addr3 = -1;
    uint64_t prev_i_Start_a_addr3 = -1;
    uint64_t prev_i_End_a_addr3 = -1;
    uint64_t prev_j_Start_a_addr3 = -1;
    uint64_t prev_j_End_a_addr3 = -1;
    uint64_t prev_cnt_a_addr4 = -1;
    uint64_t prev_i_Start_a_addr4 = -1;
    uint64_t prev_i_End_a_addr4 = -1;
    uint64_t prev_j_Start_a_addr4 = -1;
    uint64_t prev_j_End_a_addr4 = -1;
    /* Generating search reuse init code (different loops) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1025 - 1) + 1;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1025 - 1) + 1;
            j_Start = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_a_addr0 != -1) {
            if ( calAddra_addr0( i_Start - prev_i_Start_a_addr0 + prev_i_End_a_addr0, j_Start - prev_j_Start_a_addr0 + prev_j_End_a_addr0) == calAddra_addr3(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr0;
                rtHistoCal(prev_cnt_a_addr0);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr1 != -1) {
            if ( calAddra_addr1( i_Start - prev_i_Start_a_addr1 + prev_i_End_a_addr1, j_Start - prev_j_Start_a_addr1 + prev_j_End_a_addr1) == calAddra_addr3(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr1;
                rtHistoCal(prev_cnt_a_addr1);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr2 != -1) {
            if ( calAddra_addr2( i_Start - prev_i_Start_a_addr2 + prev_i_End_a_addr2, j_Start - prev_j_Start_a_addr2 + prev_j_End_a_addr2) == calAddra_addr3(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr2;
                rtHistoCal(prev_cnt_a_addr2);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr3 != -1) {
            if ( calAddra_addr3( i_Start - prev_i_Start_a_addr3 + prev_i_End_a_addr3, j_Start - prev_j_Start_a_addr3 + prev_j_End_a_addr3) == calAddra_addr3(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr3;
                rtHistoCal(prev_cnt_a_addr3);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr4 != -1) {
            if ( calAddra_addr4( i_Start - prev_i_Start_a_addr4 + prev_i_End_a_addr4, j_Start - prev_j_Start_a_addr4 + prev_j_End_a_addr4) == calAddra_addr3(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr4;
                rtHistoCal(prev_cnt_a_addr4);
                goto PROFILE;
            }
        }
        /* Generating search reuse body code (use reuse are in different loop) */
        /* Finished search reuse */
PROFILE:
sample_num++;
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
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr0 = cnt;
                        prev_i_Start_a_addr0 = i_Start;
                        prev_i_End_a_addr0 = i;
                        prev_j_Start_a_addr0 = j_Start;
                        prev_j_End_a_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr1 = cnt;
                        prev_i_Start_a_addr1 = i_Start;
                        prev_i_End_a_addr1 = i;
                        prev_j_Start_a_addr1 = j_Start;
                        prev_j_End_a_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr2 = cnt;
                        prev_i_Start_a_addr2 = i_Start;
                        prev_i_End_a_addr2 = i;
                        prev_j_Start_a_addr2 = j_Start;
                        prev_j_End_a_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr3 = cnt;
                        prev_i_Start_a_addr3 = i_Start;
                        prev_i_End_a_addr3 = i;
                        prev_j_Start_a_addr3 = j_Start;
                        prev_j_End_a_addr3 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr4 = cnt;
                        prev_i_Start_a_addr4 = i_Start;
                        prev_i_End_a_addr4 = i;
                        prev_j_Start_a_addr4 = j_Start;
                        prev_j_End_a_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        no_reuse_num++;
        if (actural != 0 && actural != pred) mis_pred_num++;
        actural = 0;
EndSample:
        s++;
        }
std::cout << "SN " << sample_num << " ";
std::cout << "NSR " << sample_num - no_reuse_num << " ";
std::cout << "RSL " << reuse_in_same_loop << " ";
std::cout << "RDL " << reuse_in_diff_loops << " ";
std::cout << "PN " << pred_num << " ";
std::cout << "PSL " << pred_sl_num << " ";
std::cout << "PDL " << pred_dl_num << " ";
std::cout << "MPN " << mis_pred_num << std::endl;
}
void ref_a_addr4() {
    /* Generating profile counter */
    uint64_t pred_num = 0;
    uint64_t pred_sl_num = 0;
    uint64_t pred_dl_num = 0;
    uint64_t sample_num = 0;
    uint64_t reuse_in_same_loop = 0;
    uint64_t reuse_in_diff_loops = 0;
    uint64_t no_reuse_num = 0;
    uint64_t mis_pred_num = 0;
    uint64_t pred = 0;
    uint64_t actural = 0;
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_a_addr0 = -1;
    uint64_t prev_i_Start_a_addr0 = -1;
    uint64_t prev_i_End_a_addr0 = -1;
    uint64_t prev_j_Start_a_addr0 = -1;
    uint64_t prev_j_End_a_addr0 = -1;
    uint64_t prev_cnt_a_addr1 = -1;
    uint64_t prev_i_Start_a_addr1 = -1;
    uint64_t prev_i_End_a_addr1 = -1;
    uint64_t prev_j_Start_a_addr1 = -1;
    uint64_t prev_j_End_a_addr1 = -1;
    uint64_t prev_cnt_a_addr2 = -1;
    uint64_t prev_i_Start_a_addr2 = -1;
    uint64_t prev_i_End_a_addr2 = -1;
    uint64_t prev_j_Start_a_addr2 = -1;
    uint64_t prev_j_End_a_addr2 = -1;
    uint64_t prev_cnt_a_addr3 = -1;
    uint64_t prev_i_Start_a_addr3 = -1;
    uint64_t prev_i_End_a_addr3 = -1;
    uint64_t prev_j_Start_a_addr3 = -1;
    uint64_t prev_j_End_a_addr3 = -1;
    uint64_t prev_cnt_a_addr4 = -1;
    uint64_t prev_i_Start_a_addr4 = -1;
    uint64_t prev_i_End_a_addr4 = -1;
    uint64_t prev_j_Start_a_addr4 = -1;
    uint64_t prev_j_End_a_addr4 = -1;
    /* Generating search reuse init code (different loops) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1025 - 1) + 1;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1025 - 1) + 1;
            j_Start = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_a_addr0 != -1) {
            if ( calAddra_addr0( i_Start - prev_i_Start_a_addr0 + prev_i_End_a_addr0, j_Start - prev_j_Start_a_addr0 + prev_j_End_a_addr0) == calAddra_addr4(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr0;
                rtHistoCal(prev_cnt_a_addr0);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr1 != -1) {
            if ( calAddra_addr1( i_Start - prev_i_Start_a_addr1 + prev_i_End_a_addr1, j_Start - prev_j_Start_a_addr1 + prev_j_End_a_addr1) == calAddra_addr4(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr1;
                rtHistoCal(prev_cnt_a_addr1);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr2 != -1) {
            if ( calAddra_addr2( i_Start - prev_i_Start_a_addr2 + prev_i_End_a_addr2, j_Start - prev_j_Start_a_addr2 + prev_j_End_a_addr2) == calAddra_addr4(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr2;
                rtHistoCal(prev_cnt_a_addr2);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr3 != -1) {
            if ( calAddra_addr3( i_Start - prev_i_Start_a_addr3 + prev_i_End_a_addr3, j_Start - prev_j_Start_a_addr3 + prev_j_End_a_addr3) == calAddra_addr4(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr3;
                rtHistoCal(prev_cnt_a_addr3);
                goto PROFILE;
            }
        }
        if ( prev_cnt_a_addr4 != -1) {
            if ( calAddra_addr4( i_Start - prev_i_Start_a_addr4 + prev_i_End_a_addr4, j_Start - prev_j_Start_a_addr4 + prev_j_End_a_addr4) == calAddra_addr4(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_a_addr4;
                rtHistoCal(prev_cnt_a_addr4);
                goto PROFILE;
            }
        }
        /* Generating search reuse body code (use reuse are in different loop) */
        /* Finished search reuse */
PROFILE:
sample_num++;
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
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr0 = cnt;
                        prev_i_Start_a_addr0 = i_Start;
                        prev_i_End_a_addr0 = i;
                        prev_j_Start_a_addr0 = j_Start;
                        prev_j_End_a_addr0 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr1( i, j) == calAddra_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr1 = cnt;
                        prev_i_Start_a_addr1 = i_Start;
                        prev_i_End_a_addr1 = i;
                        prev_j_Start_a_addr1 = j_Start;
                        prev_j_End_a_addr1 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr2( i, j) == calAddra_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr2 = cnt;
                        prev_i_Start_a_addr2 = i_Start;
                        prev_i_End_a_addr2 = i;
                        prev_j_Start_a_addr2 = j_Start;
                        prev_j_End_a_addr2 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr3( i, j) == calAddra_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr3 = cnt;
                        prev_i_Start_a_addr3 = i_Start;
                        prev_i_End_a_addr3 = i;
                        prev_j_Start_a_addr3 = j_Start;
                        prev_j_End_a_addr3 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddra_addr4( i, j) == calAddra_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_a_addr4 = cnt;
                        prev_i_Start_a_addr4 = i_Start;
                        prev_i_End_a_addr4 = i;
                        prev_j_Start_a_addr4 = j_Start;
                        prev_j_End_a_addr4 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        no_reuse_num++;
        if (actural != 0 && actural != pred) mis_pred_num++;
        actural = 0;
EndSample:
        s++;
        }
std::cout << "SN " << sample_num << " ";
std::cout << "NSR " << sample_num - no_reuse_num << " ";
std::cout << "RSL " << reuse_in_same_loop << " ";
std::cout << "RDL " << reuse_in_diff_loops << " ";
std::cout << "PN " << pred_num << " ";
std::cout << "PSL " << pred_sl_num << " ";
std::cout << "PDL " << pred_dl_num << " ";
std::cout << "MPN " << mis_pred_num << std::endl;
}
void ref_b_addr0() {
    /* Generating profile counter */
    uint64_t pred_num = 0;
    uint64_t pred_sl_num = 0;
    uint64_t pred_dl_num = 0;
    uint64_t sample_num = 0;
    uint64_t reuse_in_same_loop = 0;
    uint64_t reuse_in_diff_loops = 0;
    uint64_t no_reuse_num = 0;
    uint64_t mis_pred_num = 0;
    uint64_t pred = 0;
    uint64_t actural = 0;
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_b_addr0 = -1;
    uint64_t prev_i_Start_b_addr0 = -1;
    uint64_t prev_i_End_b_addr0 = -1;
    uint64_t prev_j_Start_b_addr0 = -1;
    uint64_t prev_j_End_b_addr0 = -1;
    /* Generating search reuse init code (different loops) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1025 - 1) + 1;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1025 - 1) + 1;
            j_Start = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_b_addr0 != -1) {
            if ( calAddrb_addr0( i_Start - prev_i_Start_b_addr0 + prev_i_End_b_addr0, j_Start - prev_j_Start_b_addr0 + prev_j_End_b_addr0) == calAddrb_addr0(i_Start, j_Start)) {
                pred_num++;
                pred_sl_num++;
                pred = prev_cnt_b_addr0;
                rtHistoCal(prev_cnt_b_addr0);
                goto PROFILE;
            }
        }
        /* Generating search reuse body code (use reuse are in different loop) */
        /* Finished search reuse */
PROFILE:
sample_num++;
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
                        actural = cnt;
                        reuse_in_same_loop++;
                        prev_cnt_b_addr0 = cnt;
                        prev_i_Start_b_addr0 = i_Start;
                        prev_i_End_b_addr0 = i;
                        prev_j_Start_b_addr0 = j_Start;
                        prev_j_End_b_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
        }
        }
        no_reuse_num++;
        if (actural != 0 && actural != pred) mis_pred_num++;
        actural = 0;
EndSample:
        s++;
        }
std::cout << "SN " << sample_num << " ";
std::cout << "NSR " << sample_num - no_reuse_num << " ";
std::cout << "RSL " << reuse_in_same_loop << " ";
std::cout << "RDL " << reuse_in_diff_loops << " ";
std::cout << "PN " << pred_num << " ";
std::cout << "PSL " << pred_sl_num << " ";
std::cout << "PDL " << pred_dl_num << " ";
std::cout << "MPN " << mis_pred_num << std::endl;
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
 /* Start to analyze function:  
stencil */ 
