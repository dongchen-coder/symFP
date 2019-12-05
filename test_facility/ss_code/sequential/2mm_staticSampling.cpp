
 /* Start to analysis array index
Array index info: Total number of references: 11
tmp.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
C.addr ((k * 1024) + j)
D.addr ((i * 1024) + j)
D.addr ((i * 1024) + j)
B.addr ((k * 1024) + j)
tmp.addr ((i * 1024) + j)
tmp.addr ((i * 1024) + j)
D.addr ((i * 1024) + j)
D.addr ((i * 1024) + j)
tmp.addr ((i * 1024) + k)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %tmp
double* %A
double* %B
double* %C
double* %D
double %alpha
double %beta

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access tmp.addr ((i * 1024) + j)
------k
------Loop Bound: (0, 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 1024) + k)
--------array access B.addr ((k * 1024) + j)
--------array access tmp.addr ((i * 1024) + j)
--------array access tmp.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access D.addr ((i * 1024) + j)
------array access D.addr ((i * 1024) + j)
------k
------Loop Bound: (0, 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access tmp.addr ((i * 1024) + k)
--------array access C.addr ((k * 1024) + j)
--------array access D.addr ((i * 1024) + j)
--------array access D.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 51
------Sample number: 2621
--------Sample number: 134217
----Sample number: 51
------Sample number: 2621
--------Sample number: 134217
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
#include <cmath>
#ifdef PAPI_TIMER
#  include <chrono>
#endif
using namespace std;
#ifdef PAPI_TIMER
using namespace std::chrono;
#endif
map<uint64_t, map<uint64_t, uint64_t>* > RI;
map<uint64_t, map<uint64_t, double>* > hits;
map<uint64_t, map<uint64_t, double>* > costs;
map<uint64_t, double> sampledCnt;
map<uint64_t, double> accessRatio;
map<uint64_t, uint64_t> Lease;
void rtHistoCal(uint64_t ri, uint64_t ref_id) {
    if (RI.find(ref_id) != RI.end()) {
        if ((*RI[ref_id]).find(ri) != (*RI[ref_id]).end()) {
            (*RI[ref_id])[ri] ++;
        } else {
            (*RI[ref_id])[ri] = 1;
        }
    } else {
        RI[ref_id] = new map<uint64_t, uint64_t>;
        (*RI[ref_id])[ri] = 1;
    }

    // Init leases to all references to be 0
    if (Lease.find(ref_id) == Lease.end()) {
        Lease[ref_id] = 0;
    }
    return;
}
void accessRatioCal() {
    double total_access_cnt = 0;

    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
        for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            total_access_cnt += ri_it->second;
        }
    }
    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
        double ref_access_cnt = 0;
        for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
           ref_access_cnt += ri_it->second;
       }
        sampledCnt[ref_it->first] = ref_access_cnt;
        accessRatio[ref_it->first] = ref_access_cnt / total_access_cnt;
    }
}
void initHitsCosts() {
    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
        hits[ref_it->first] = new map<uint64_t, double>;
        costs[ref_it->first] = new map<uint64_t, double>;
        (*hits[ref_it->first])[0] = 0;
        uint64_t total_hits = 0;
        (*costs[ref_it->first])[0] = 0;
        uint64_t total_cnt = 0;
        for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            total_cnt += ri_it->second;
        }
        uint64_t pre_lease = 0;
        uint64_t pre_cost = 0;
        for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            total_hits += ri_it->second;
            (*hits[ref_it->first])[ri_it->first] = total_hits;
            (*costs[ref_it->first])[ri_it->first] =  pre_cost + (ri_it->first - pre_lease) * total_cnt;
            total_cnt -= ri_it->second;
            pre_cost = (*costs[ref_it->first])[ri_it->first];
            pre_lease = ri_it->first;
        }
    }
}
double getPPUC(uint64_t ref_id, uint64_t oldLease, uint64_t newLease) {
    if (hits.find(ref_id) == hits.end() || costs.find(ref_id) == costs.end()) {
        cout << "No such ref for hits/costs" << endl;
        return -1;
    }
    if (hits[ref_id]->find(newLease) == hits[ref_id]->end() || costs[ref_id]->find(newLease) == costs[ref_id]->end()) {
        cout << "No RI/Newlease " << newLease << " for ref " << ref_id << endl;
        return -1;
    }
    if (hits[ref_id]->find(oldLease) == hits[ref_id]->end() || costs[ref_id]->find(oldLease) == costs[ref_id]->end()) {
        if (hits[ref_id]->find(oldLease) == hits[ref_id]->end()) {
            cout << "No hits for Oldlease " << oldLease << " for ref " << ref_id << endl;
        }
        if (costs[ref_id]->find(oldLease) == costs[ref_id]->end()) {
            cout << "No costs for Oldlease " << oldLease << " for ref " << ref_id << endl;
        }
        return -1;
    }
    return double((*hits[ref_id])[newLease] - (*hits[ref_id])[oldLease]) / ((*costs[ref_id])[newLease] - (*costs[ref_id])[oldLease]);
}
void getMaxPPUC(bool*finished, uint64_t* ref_to_assign, uint64_t* newLease) {
    double maxPPUC = -1;
    uint64_t bestRef = -1;
    uint64_t bestLease = -1;
    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
        for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            if (ri_it->first > Lease[ref_it->first]) {
                double ppuc = getPPUC(ref_it->first, Lease[ref_it->first], ri_it->first);
                if (ppuc > maxPPUC) {
                    maxPPUC = ppuc;
                    bestRef = ref_it->first;
                    bestLease = ri_it->first;
                }
            }
        }
    }
    if (maxPPUC != -1) {
        *finished = false;
        *ref_to_assign = bestRef;
        *newLease = bestLease;
    } else {
        *finished = true;
    }
    return;
}
void dumpRI() {
    uint64_t total_number_of_ri = 0;
    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
        std::set<uint64_t> riset;
        for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            cout << "Ref " << ref_it->first << " RI " << ri_it->first << " CNT " << ri_it->second << endl;
            riset.insert(ri_it->first);
        }
        cout << "Ref " << ref_it->first << " RISETSIZE " << riset.size() << endl;
        total_number_of_ri += riset.size();
    }
    cout << "Average RISETSIZE for each reference " << double(total_number_of_ri) / RI.size() << endl;
}
void RL_main(uint64_t CacheSize) {
    initHitsCosts();
    accessRatioCal();
    double totalCost = 0;
    double totalHitRatio = 0;
    double targetCost = CacheSize;
#ifdef DEBUG
    dumpRI();
#endif
    while(true) {
        bool finished = false;
        uint64_t ref_to_assign;
        uint64_t newLease;
        getMaxPPUC(&finished, &ref_to_assign, &newLease);
        if (finished == false) {
            totalCost += ((*costs[ref_to_assign])[newLease] - (*costs[ref_to_assign])[Lease[ref_to_assign]]) / sampledCnt[ref_to_assign] * accessRatio[ref_to_assign];
            totalHitRatio += ((*hits[ref_to_assign])[newLease] - (*hits[ref_to_assign])[Lease[ref_to_assign]]) / sampledCnt[ref_to_assign] * accessRatio[ref_to_assign];
            Lease[ref_to_assign] = newLease;
            cout << "Assign lease " << newLease << " to ref " << ref_to_assign << " avg cache size " << totalCost  << " miss ratio " << 1 - totalHitRatio << endl;
        } else {
            break;
        }
        if (totalCost < targetCost && targetCost != 0) {
            break;
        }
    }
    return;
}
/* tmp_addr ((i * 1024) + j) 0 */
int calAddrtmp_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + k) 1 */
int calAddrA_addr1( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* B_addr ((k * 1024) + j) 2 */
int calAddrB_addr2( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
/* tmp_addr ((i * 1024) + j) 3 */
int calAddrtmp_addr3( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* tmp_addr ((i * 1024) + j) 4 */
int calAddrtmp_addr4( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 1024) + j) 5 */
int calAddrD_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 1024) + j) 6 */
int calAddrD_addr6( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* tmp_addr ((i * 1024) + k) 7 */
int calAddrtmp_addr7( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* C_addr ((k * 1024) + j) 8 */
int calAddrC_addr8( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 1024) + j) 9 */
int calAddrD_addr9( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 1024) + j) 10 */
int calAddrD_addr10( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_tmp_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_tmp_addr0 = -1;
    uint64_t prev_i_Start_tmp_addr0 = -1;
    uint64_t prev_i_End_tmp_addr0 = -1;
    uint64_t prev_j_Start_tmp_addr0 = -1;
    uint64_t prev_j_End_tmp_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2621;) {
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
        if ( prev_cnt_tmp_addr0 != -1) {
            if ( calAddrtmp_addr0( i_Start - prev_i_Start_tmp_addr0 + prev_i_End_tmp_addr0, j_Start - prev_j_Start_tmp_addr0 + prev_j_End_tmp_addr0) == calAddrtmp_addr0(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_tmp_addr0, 0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr0( i, j) == calAddrtmp_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt, 0);
                        prev_cnt_tmp_addr0 = cnt;
                        prev_i_Start_tmp_addr0 = i_Start;
                        prev_i_End_tmp_addr0 = i;
                        prev_j_Start_tmp_addr0 = j_Start;
                        prev_j_End_tmp_addr0 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr3( i, j, k) == calAddrtmp_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt, 0);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr4( i, j, k) == calAddrtmp_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt, 0);
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
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr7( i, j, k) == calAddrtmp_addr0(i_Start, j_Start)) {
                            rtHistoCal(cnt, 0);
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
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_k_Start_A_addr1 = -1;
    uint64_t prev_k_End_A_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 134217;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1, k_Start - prev_k_Start_A_addr1 + prev_k_End_A_addr1) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr1, 1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 1);
                            prev_cnt_A_addr1 = cnt;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            prev_k_Start_A_addr1 = k_Start;
                            prev_k_End_A_addr1 = k;
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
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k++) {
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
void ref_C_addr8() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_C_addr8 = -1;
    uint64_t prev_i_Start_C_addr8 = -1;
    uint64_t prev_i_End_C_addr8 = -1;
    uint64_t prev_j_Start_C_addr8 = -1;
    uint64_t prev_j_End_C_addr8 = -1;
    uint64_t prev_k_Start_C_addr8 = -1;
    uint64_t prev_k_End_C_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 134217;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_C_addr8 != -1) {
            if ( calAddrC_addr8( i_Start - prev_i_Start_C_addr8 + prev_i_End_C_addr8, j_Start - prev_j_Start_C_addr8 + prev_j_End_C_addr8, k_Start - prev_k_Start_C_addr8 + prev_k_End_C_addr8) == calAddrC_addr8(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_C_addr8, 8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrC_addr8( i, j, k) == calAddrC_addr8(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 8);
                            prev_cnt_C_addr8 = cnt;
                            prev_i_Start_C_addr8 = i_Start;
                            prev_i_End_C_addr8 = i;
                            prev_j_Start_C_addr8 = j_Start;
                            prev_j_End_C_addr8 = j;
                            prev_k_Start_C_addr8 = k_Start;
                            prev_k_End_C_addr8 = k;
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
        }
        }
EndSample:
        s++;
        }
}
void ref_D_addr9() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_D_addr9 = -1;
    uint64_t prev_i_Start_D_addr9 = -1;
    uint64_t prev_i_End_D_addr9 = -1;
    uint64_t prev_j_Start_D_addr9 = -1;
    uint64_t prev_j_End_D_addr9 = -1;
    uint64_t prev_k_Start_D_addr9 = -1;
    uint64_t prev_k_End_D_addr9 = -1;
    uint64_t prev_cnt_D_addr10 = -1;
    uint64_t prev_i_Start_D_addr10 = -1;
    uint64_t prev_i_End_D_addr10 = -1;
    uint64_t prev_j_Start_D_addr10 = -1;
    uint64_t prev_j_End_D_addr10 = -1;
    uint64_t prev_k_Start_D_addr10 = -1;
    uint64_t prev_k_End_D_addr10 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 134217;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_D_addr9 != -1) {
            if ( calAddrD_addr9( i_Start - prev_i_Start_D_addr9 + prev_i_End_D_addr9, j_Start - prev_j_Start_D_addr9 + prev_j_End_D_addr9, k_Start - prev_k_Start_D_addr9 + prev_k_End_D_addr9) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_D_addr9, 9);
                goto EndSample;
            }
        }
        if ( prev_cnt_D_addr10 != -1) {
            if ( calAddrD_addr10( i_Start - prev_i_Start_D_addr10 + prev_i_End_D_addr10, j_Start - prev_j_Start_D_addr10 + prev_j_End_D_addr10, k_Start - prev_k_Start_D_addr10 + prev_k_End_D_addr10) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_D_addr10, 10);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr5( i, j) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt, 9);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr6( i, j) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt, 9);
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr9( i, j, k) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 9);
                            prev_cnt_D_addr9 = cnt;
                            prev_i_Start_D_addr9 = i_Start;
                            prev_i_End_D_addr9 = i;
                            prev_j_Start_D_addr9 = j_Start;
                            prev_j_End_D_addr9 = j;
                            prev_k_Start_D_addr9 = k_Start;
                            prev_k_End_D_addr9 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr10( i, j, k) == calAddrD_addr9(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 9);
                            prev_cnt_D_addr10 = cnt;
                            prev_i_Start_D_addr10 = i_Start;
                            prev_i_End_D_addr10 = i;
                            prev_j_Start_D_addr10 = j_Start;
                            prev_j_End_D_addr10 = j;
                            prev_k_Start_D_addr10 = k_Start;
                            prev_k_End_D_addr10 = k;
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
void ref_D_addr10() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_D_addr9 = -1;
    uint64_t prev_i_Start_D_addr9 = -1;
    uint64_t prev_i_End_D_addr9 = -1;
    uint64_t prev_j_Start_D_addr9 = -1;
    uint64_t prev_j_End_D_addr9 = -1;
    uint64_t prev_k_Start_D_addr9 = -1;
    uint64_t prev_k_End_D_addr9 = -1;
    uint64_t prev_cnt_D_addr10 = -1;
    uint64_t prev_i_Start_D_addr10 = -1;
    uint64_t prev_i_End_D_addr10 = -1;
    uint64_t prev_j_Start_D_addr10 = -1;
    uint64_t prev_j_End_D_addr10 = -1;
    uint64_t prev_k_Start_D_addr10 = -1;
    uint64_t prev_k_End_D_addr10 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 134217;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_D_addr9 != -1) {
            if ( calAddrD_addr9( i_Start - prev_i_Start_D_addr9 + prev_i_End_D_addr9, j_Start - prev_j_Start_D_addr9 + prev_j_End_D_addr9, k_Start - prev_k_Start_D_addr9 + prev_k_End_D_addr9) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_D_addr9, 9);
                goto EndSample;
            }
        }
        if ( prev_cnt_D_addr10 != -1) {
            if ( calAddrD_addr10( i_Start - prev_i_Start_D_addr10 + prev_i_End_D_addr10, j_Start - prev_j_Start_D_addr10 + prev_j_End_D_addr10, k_Start - prev_k_Start_D_addr10 + prev_k_End_D_addr10) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_D_addr10, 10);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr5( i, j) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt, 10);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr6( i, j) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt, 10);
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr9( i, j, k) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 10);
                            prev_cnt_D_addr9 = cnt;
                            prev_i_Start_D_addr9 = i_Start;
                            prev_i_End_D_addr9 = i;
                            prev_j_Start_D_addr9 = j_Start;
                            prev_j_End_D_addr9 = j;
                            prev_k_Start_D_addr9 = k_Start;
                            prev_k_End_D_addr9 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr10( i, j, k) == calAddrD_addr10(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 10);
                            prev_cnt_D_addr10 = cnt;
                            prev_i_Start_D_addr10 = i_Start;
                            prev_i_End_D_addr10 = i;
                            prev_j_Start_D_addr10 = j_Start;
                            prev_j_End_D_addr10 = j;
                            prev_k_Start_D_addr10 = k_Start;
                            prev_k_End_D_addr10 = k;
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
void ref_B_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_B_addr2 = -1;
    uint64_t prev_i_Start_B_addr2 = -1;
    uint64_t prev_i_End_B_addr2 = -1;
    uint64_t prev_j_Start_B_addr2 = -1;
    uint64_t prev_j_End_B_addr2 = -1;
    uint64_t prev_k_Start_B_addr2 = -1;
    uint64_t prev_k_End_B_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 134217;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_B_addr2 != -1) {
            if ( calAddrB_addr2( i_Start - prev_i_Start_B_addr2 + prev_i_End_B_addr2, j_Start - prev_j_Start_B_addr2 + prev_j_End_B_addr2, k_Start - prev_k_Start_B_addr2 + prev_k_End_B_addr2) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_B_addr2, 2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrB_addr2( i, j, k) == calAddrB_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 2);
                            prev_cnt_B_addr2 = cnt;
                            prev_i_Start_B_addr2 = i_Start;
                            prev_i_End_B_addr2 = i;
                            prev_j_Start_B_addr2 = j_Start;
                            prev_j_End_B_addr2 = j;
                            prev_k_Start_B_addr2 = k_Start;
                            prev_k_End_B_addr2 = k;
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
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k++) {
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
void ref_tmp_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_tmp_addr3 = -1;
    uint64_t prev_i_Start_tmp_addr3 = -1;
    uint64_t prev_i_End_tmp_addr3 = -1;
    uint64_t prev_j_Start_tmp_addr3 = -1;
    uint64_t prev_j_End_tmp_addr3 = -1;
    uint64_t prev_k_Start_tmp_addr3 = -1;
    uint64_t prev_k_End_tmp_addr3 = -1;
    uint64_t prev_cnt_tmp_addr4 = -1;
    uint64_t prev_i_Start_tmp_addr4 = -1;
    uint64_t prev_i_End_tmp_addr4 = -1;
    uint64_t prev_j_Start_tmp_addr4 = -1;
    uint64_t prev_j_End_tmp_addr4 = -1;
    uint64_t prev_k_Start_tmp_addr4 = -1;
    uint64_t prev_k_End_tmp_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 134217;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_tmp_addr3 != -1) {
            if ( calAddrtmp_addr3( i_Start - prev_i_Start_tmp_addr3 + prev_i_End_tmp_addr3, j_Start - prev_j_Start_tmp_addr3 + prev_j_End_tmp_addr3, k_Start - prev_k_Start_tmp_addr3 + prev_k_End_tmp_addr3) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_tmp_addr3, 3);
                goto EndSample;
            }
        }
        if ( prev_cnt_tmp_addr4 != -1) {
            if ( calAddrtmp_addr4( i_Start - prev_i_Start_tmp_addr4 + prev_i_End_tmp_addr4, j_Start - prev_j_Start_tmp_addr4 + prev_j_End_tmp_addr4, k_Start - prev_k_Start_tmp_addr4 + prev_k_End_tmp_addr4) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_tmp_addr4, 4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr0( i, j) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt, 3);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr3( i, j, k) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 3);
                            prev_cnt_tmp_addr3 = cnt;
                            prev_i_Start_tmp_addr3 = i_Start;
                            prev_i_End_tmp_addr3 = i;
                            prev_j_Start_tmp_addr3 = j_Start;
                            prev_j_End_tmp_addr3 = j;
                            prev_k_Start_tmp_addr3 = k_Start;
                            prev_k_End_tmp_addr3 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr4( i, j, k) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 3);
                            prev_cnt_tmp_addr4 = cnt;
                            prev_i_Start_tmp_addr4 = i_Start;
                            prev_i_End_tmp_addr4 = i;
                            prev_j_Start_tmp_addr4 = j_Start;
                            prev_j_End_tmp_addr4 = j;
                            prev_k_Start_tmp_addr4 = k_Start;
                            prev_k_End_tmp_addr4 = k;
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
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr7( i, j, k) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 3);
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
        }
        }
EndSample:
        s++;
        }
}
void ref_tmp_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_tmp_addr3 = -1;
    uint64_t prev_i_Start_tmp_addr3 = -1;
    uint64_t prev_i_End_tmp_addr3 = -1;
    uint64_t prev_j_Start_tmp_addr3 = -1;
    uint64_t prev_j_End_tmp_addr3 = -1;
    uint64_t prev_k_Start_tmp_addr3 = -1;
    uint64_t prev_k_End_tmp_addr3 = -1;
    uint64_t prev_cnt_tmp_addr4 = -1;
    uint64_t prev_i_Start_tmp_addr4 = -1;
    uint64_t prev_i_End_tmp_addr4 = -1;
    uint64_t prev_j_Start_tmp_addr4 = -1;
    uint64_t prev_j_End_tmp_addr4 = -1;
    uint64_t prev_k_Start_tmp_addr4 = -1;
    uint64_t prev_k_End_tmp_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 134217;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_tmp_addr3 != -1) {
            if ( calAddrtmp_addr3( i_Start - prev_i_Start_tmp_addr3 + prev_i_End_tmp_addr3, j_Start - prev_j_Start_tmp_addr3 + prev_j_End_tmp_addr3, k_Start - prev_k_Start_tmp_addr3 + prev_k_End_tmp_addr3) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_tmp_addr3, 3);
                goto EndSample;
            }
        }
        if ( prev_cnt_tmp_addr4 != -1) {
            if ( calAddrtmp_addr4( i_Start - prev_i_Start_tmp_addr4 + prev_i_End_tmp_addr4, j_Start - prev_j_Start_tmp_addr4 + prev_j_End_tmp_addr4, k_Start - prev_k_Start_tmp_addr4 + prev_k_End_tmp_addr4) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_tmp_addr4, 4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrtmp_addr0( i, j) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt, 4);
                        goto EndSample;
                    }
                }
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr3( i, j, k) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 4);
                            prev_cnt_tmp_addr3 = cnt;
                            prev_i_Start_tmp_addr3 = i_Start;
                            prev_i_End_tmp_addr3 = i;
                            prev_j_Start_tmp_addr3 = j_Start;
                            prev_j_End_tmp_addr3 = j;
                            prev_k_Start_tmp_addr3 = k_Start;
                            prev_k_End_tmp_addr3 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr4( i, j, k) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 4);
                            prev_cnt_tmp_addr4 = cnt;
                            prev_i_Start_tmp_addr4 = i_Start;
                            prev_i_End_tmp_addr4 = i;
                            prev_j_Start_tmp_addr4 = j_Start;
                            prev_j_End_tmp_addr4 = j;
                            prev_k_Start_tmp_addr4 = k_Start;
                            prev_k_End_tmp_addr4 = k;
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
        int iLB3 = 0;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr7( i, j, k) == calAddrtmp_addr4(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 4);
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
        }
        }
EndSample:
        s++;
        }
}
void ref_D_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_D_addr5 = -1;
    uint64_t prev_i_Start_D_addr5 = -1;
    uint64_t prev_i_End_D_addr5 = -1;
    uint64_t prev_j_Start_D_addr5 = -1;
    uint64_t prev_j_End_D_addr5 = -1;
    uint64_t prev_cnt_D_addr6 = -1;
    uint64_t prev_i_Start_D_addr6 = -1;
    uint64_t prev_i_End_D_addr6 = -1;
    uint64_t prev_j_Start_D_addr6 = -1;
    uint64_t prev_j_End_D_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2621;) {
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
        if ( prev_cnt_D_addr5 != -1) {
            if ( calAddrD_addr5( i_Start - prev_i_Start_D_addr5 + prev_i_End_D_addr5, j_Start - prev_j_Start_D_addr5 + prev_j_End_D_addr5) == calAddrD_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_D_addr5, 5);
                goto EndSample;
            }
        }
        if ( prev_cnt_D_addr6 != -1) {
            if ( calAddrD_addr6( i_Start - prev_i_Start_D_addr6 + prev_i_End_D_addr6, j_Start - prev_j_Start_D_addr6 + prev_j_End_D_addr6) == calAddrD_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_D_addr6, 6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr5( i, j) == calAddrD_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt, 5);
                        prev_cnt_D_addr5 = cnt;
                        prev_i_Start_D_addr5 = i_Start;
                        prev_i_End_D_addr5 = i;
                        prev_j_Start_D_addr5 = j_Start;
                        prev_j_End_D_addr5 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr6( i, j) == calAddrD_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt, 5);
                        prev_cnt_D_addr6 = cnt;
                        prev_i_Start_D_addr6 = i_Start;
                        prev_i_End_D_addr6 = i;
                        prev_j_Start_D_addr6 = j_Start;
                        prev_j_End_D_addr6 = j;
                        goto EndSample;
                    }
                }
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr9( i, j, k) == calAddrD_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr10( i, j, k) == calAddrD_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
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
void ref_D_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_D_addr5 = -1;
    uint64_t prev_i_Start_D_addr5 = -1;
    uint64_t prev_i_End_D_addr5 = -1;
    uint64_t prev_j_Start_D_addr5 = -1;
    uint64_t prev_j_End_D_addr5 = -1;
    uint64_t prev_cnt_D_addr6 = -1;
    uint64_t prev_i_Start_D_addr6 = -1;
    uint64_t prev_i_End_D_addr6 = -1;
    uint64_t prev_j_Start_D_addr6 = -1;
    uint64_t prev_j_End_D_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2621;) {
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
        if ( prev_cnt_D_addr5 != -1) {
            if ( calAddrD_addr5( i_Start - prev_i_Start_D_addr5 + prev_i_End_D_addr5, j_Start - prev_j_Start_D_addr5 + prev_j_End_D_addr5) == calAddrD_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_D_addr5, 5);
                goto EndSample;
            }
        }
        if ( prev_cnt_D_addr6 != -1) {
            if ( calAddrD_addr6( i_Start - prev_i_Start_D_addr6 + prev_i_End_D_addr6, j_Start - prev_j_Start_D_addr6 + prev_j_End_D_addr6) == calAddrD_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_D_addr6, 6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr5( i, j) == calAddrD_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt, 6);
                        prev_cnt_D_addr5 = cnt;
                        prev_i_Start_D_addr5 = i_Start;
                        prev_i_End_D_addr5 = i;
                        prev_j_Start_D_addr5 = j_Start;
                        prev_j_End_D_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrD_addr6( i, j) == calAddrD_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt, 6);
                        prev_cnt_D_addr6 = cnt;
                        prev_i_Start_D_addr6 = i_Start;
                        prev_i_End_D_addr6 = i;
                        prev_j_Start_D_addr6 = j_Start;
                        prev_j_End_D_addr6 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int kLB5 = 0;
                for ( int k = kLB5; k < 1024; k++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr9( i, j, k) == calAddrD_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt, 6);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrD_addr10( i, j, k) == calAddrD_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt, 6);
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
void ref_tmp_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_tmp_addr7 = -1;
    uint64_t prev_i_Start_tmp_addr7 = -1;
    uint64_t prev_i_End_tmp_addr7 = -1;
    uint64_t prev_j_Start_tmp_addr7 = -1;
    uint64_t prev_j_End_tmp_addr7 = -1;
    uint64_t prev_k_Start_tmp_addr7 = -1;
    uint64_t prev_k_End_tmp_addr7 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 134217;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_tmp_addr7 != -1) {
            if ( calAddrtmp_addr7( i_Start - prev_i_Start_tmp_addr7 + prev_i_End_tmp_addr7, j_Start - prev_j_Start_tmp_addr7 + prev_j_End_tmp_addr7, k_Start - prev_k_Start_tmp_addr7 + prev_k_End_tmp_addr7) == calAddrtmp_addr7(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_tmp_addr7, 7);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 1024; i++) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int kLB5 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB5 = k_Start;
                }
                for ( int k = kLB5; k < 1024; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrtmp_addr7( i, j, k) == calAddrtmp_addr7(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt, 7);
                            prev_cnt_tmp_addr7 = cnt;
                            prev_i_Start_tmp_addr7 = i_Start;
                            prev_i_End_tmp_addr7 = i;
                            prev_j_Start_tmp_addr7 = j_Start;
                            prev_j_End_tmp_addr7 = j;
                            prev_k_Start_tmp_addr7 = k_Start;
                            prev_k_End_tmp_addr7 = k;
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
        }
        }
EndSample:
        s++;
        }
}
int main() {
#ifdef PAPI_TIMER
    // Get starting timepoint
    auto start = high_resolution_clock::now();
#endif
    ref_tmp_addr0();
    ref_A_addr1();
    ref_C_addr8();
    ref_D_addr9();
    ref_D_addr10();
    ref_B_addr2();
    ref_tmp_addr3();
    ref_tmp_addr4();
    ref_D_addr5();
    ref_D_addr6();
    ref_tmp_addr7();
#ifdef PAPI_TIMER
// Get ending timepoint
    auto stop = high_resolution_clock::now(); 
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
     cout << "Time taken by SPS:  " << duration.count() << endl; 
#endif
#ifdef PAPI_TIMER
    // Get starting timepoint
    start = high_resolution_clock::now();
#endif
    RL_main(0);
#ifdef PAPI_TIMER
// Get ending timepoint
    stop = high_resolution_clock::now(); 
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    duration = duration_cast<microseconds>(stop - start);
     cout << "Time taken by CARL:  " << duration.count() << endl; 
#endif
    return 0;
}
 /* Analyze function: mm2 */ 
