
 /* Start to analysis array index
Array index info: Total number of references: 3
pC.addr ((32 * i) + j)
pA.addr ((32 * i) + k)
pB.addr ((32 * k) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %pA
double* %pB
double* %pC

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--ii
--Loop Bound: (0, 32)
--Loop inc: (ii + 8)
--Loop predicate: <
----jj
----Loop Bound: (0, 32)
----Loop inc: (jj + 8)
----Loop predicate: <
------kk
------Loop Bound: (0, 32)
------Loop inc: (kk + 8)
------Loop predicate: <
--------i
--------Loop Bound: (ii, (ii + 8))
--------Loop inc: (i + 1)
--------Loop predicate: <
----------j
----------Loop Bound: (jj, (jj + 8))
----------Loop inc: (j + 1)
----------Loop predicate: <
------------k
------------Loop Bound: (kk, (kk + 8))
------------Loop inc: (k + 1)
------------Loop predicate: <
--------------array access pA.addr ((32 * i) + k)
--------------array access pB.addr ((32 * k) + j)
------------array access pC.addr ((32 * i) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 0 0 0 0 
Dump stride: 8 8 8 1 
init counter: 0 0 0 0 0 
Dump stride: 8 8 8 1 1 
init counter: 0 0 0 0 0 0 
Dump stride: 8 8 8 1 1 1 
Dump tree:
----Sample number: 4
------Sample number: 16
--------Sample number: 64
----------Sample number: 512
------------Sample number: 4096
--------------Sample number: 32768
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
using namespace  chrono;
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
         set<uint64_t> riset;
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
/* pA_addr ((32 * i) + k) 0 */
int calAddrpA_addr0( int ii, int jj, int kk, int i, int j, int k) {
    int result = (((32 * i) + k)) * 8 / 64;
    return result;
}
/* pB_addr ((32 * k) + j) 1 */
int calAddrpB_addr1( int ii, int jj, int kk, int i, int j, int k) {
    int result = (((32 * k) + j)) * 8 / 64;
    return result;
}
/* pC_addr ((32 * i) + j) 2 */
int calAddrpC_addr2( int ii, int jj, int kk, int i, int j) {
    int result = (((32 * i) + j)) * 8 / 64;
    return result;
}
void ref_pC_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 4096;) {
SAMPLE:
        int ii_Start = rand() % (32 - 0) + 0;
        if ( (32 - 0) == 0) goto SAMPLE;
        int jj_Start = rand() % (32 - 0) + 0;
        if ( (32 - 0) == 0) goto SAMPLE;
        int kk_Start = rand() % (32 - 0) + 0;
        if ( ((ii_Start + 8) - ii_Start) == 0) goto SAMPLE;
        int i_Start = rand() % ((ii_Start + 8) - ii_Start) + ii_Start;
        if ( ((jj_Start + 8) - jj_Start) == 0) goto SAMPLE;
        int j_Start = rand() % ((jj_Start + 8) - jj_Start) + jj_Start;
        string idx_string =  to_string(ii_Start) + "_" +  to_string(jj_Start) + "_" +  to_string(kk_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iiLB0 = ii_Start;
        for ( int ii = iiLB0; ii < 32; ii++) {
            {
            int jjLB1 = 0;
            if ( ii == ii_Start ) {
                jjLB1 = jj_Start;
            }
            for ( int jj = jjLB1; jj < 32; jj++) {
                {
                int kkLB2 = 0;
                if ( ii == ii_Start && jj == jj_Start ) {
                    kkLB2 = kk_Start;
                }
                for ( int kk = kkLB2; kk < 32; kk++) {
                    {
                    int iLB3 = ii;
                    if ( ii == ii_Start && jj == jj_Start && kk == kk_Start ) {
                        iLB3 = i_Start;
                    }
                    for ( int i = iLB3; i < (ii + 8); i++) {
                        {
                        int jLB4 = jj;
                        if ( ii == ii_Start && jj == jj_Start && kk == kk_Start && i == i_Start ) {
                            jLB4 = j_Start;
                        }
                        for ( int j = jLB4; j < (jj + 8); j++) {
                            {
                            int kLB5 = kk;
                            for ( int k = kLB5; k < (kk + 8); k++) {
                                if (cntStart == true) cnt++;
                                if (cntStart == true) cnt++;
                            }
                            }
                            if (cntStart == true) {
                                cnt++;
                                if ( calAddrpC_addr2( ii, jj, kk, i, j) == calAddrpC_addr2(ii_Start, jj_Start, kk_Start, i_Start, j_Start)) {
                                    rtHistoCal(cnt, 2);
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
        }
        }
EndSample:
        s++;
        }
}
void ref_pA_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 32768;) {
SAMPLE:
        int ii_Start = rand() % (32 - 0) + 0;
        if ( (32 - 0) == 0) goto SAMPLE;
        int jj_Start = rand() % (32 - 0) + 0;
        if ( (32 - 0) == 0) goto SAMPLE;
        int kk_Start = rand() % (32 - 0) + 0;
        if ( ((ii_Start + 8) - ii_Start) == 0) goto SAMPLE;
        int i_Start = rand() % ((ii_Start + 8) - ii_Start) + ii_Start;
        if ( ((jj_Start + 8) - jj_Start) == 0) goto SAMPLE;
        int j_Start = rand() % ((jj_Start + 8) - jj_Start) + jj_Start;
        if ( ((kk_Start + 8) - kk_Start) == 0) goto SAMPLE;
        int k_Start = rand() % ((kk_Start + 8) - kk_Start) + kk_Start;
        string idx_string =  to_string(ii_Start) + "_" +  to_string(jj_Start) + "_" +  to_string(kk_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iiLB0 = ii_Start;
        for ( int ii = iiLB0; ii < 32; ii++) {
            {
            int jjLB1 = 0;
            if ( ii == ii_Start ) {
                jjLB1 = jj_Start;
            }
            for ( int jj = jjLB1; jj < 32; jj++) {
                {
                int kkLB2 = 0;
                if ( ii == ii_Start && jj == jj_Start ) {
                    kkLB2 = kk_Start;
                }
                for ( int kk = kkLB2; kk < 32; kk++) {
                    {
                    int iLB3 = ii;
                    if ( ii == ii_Start && jj == jj_Start && kk == kk_Start ) {
                        iLB3 = i_Start;
                    }
                    for ( int i = iLB3; i < (ii + 8); i++) {
                        {
                        int jLB4 = jj;
                        if ( ii == ii_Start && jj == jj_Start && kk == kk_Start && i == i_Start ) {
                            jLB4 = j_Start;
                        }
                        for ( int j = jLB4; j < (jj + 8); j++) {
                            {
                            int kLB5 = kk;
                            if ( ii == ii_Start && jj == jj_Start && kk == kk_Start && i == i_Start && j == j_Start ) {
                                kLB5 = k_Start;
                            }
                            for ( int k = kLB5; k < (kk + 8); k++) {
                                if (cntStart == true) {
                                    cnt++;
                                    if ( calAddrpA_addr0( ii, jj, kk, i, j, k) == calAddrpA_addr0(ii_Start, jj_Start, kk_Start, i_Start, j_Start, k_Start)) {
                                        rtHistoCal(cnt, 0);
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
void ref_pB_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 32768;) {
SAMPLE:
        int ii_Start = rand() % (32 - 0) + 0;
        if ( (32 - 0) == 0) goto SAMPLE;
        int jj_Start = rand() % (32 - 0) + 0;
        if ( (32 - 0) == 0) goto SAMPLE;
        int kk_Start = rand() % (32 - 0) + 0;
        if ( ((ii_Start + 8) - ii_Start) == 0) goto SAMPLE;
        int i_Start = rand() % ((ii_Start + 8) - ii_Start) + ii_Start;
        if ( ((jj_Start + 8) - jj_Start) == 0) goto SAMPLE;
        int j_Start = rand() % ((jj_Start + 8) - jj_Start) + jj_Start;
        if ( ((kk_Start + 8) - kk_Start) == 0) goto SAMPLE;
        int k_Start = rand() % ((kk_Start + 8) - kk_Start) + kk_Start;
        string idx_string =  to_string(ii_Start) + "_" +  to_string(jj_Start) + "_" +  to_string(kk_Start) + "_" +  to_string(i_Start) + "_" +  to_string(j_Start) + "_" +  to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iiLB0 = ii_Start;
        for ( int ii = iiLB0; ii < 32; ii++) {
            {
            int jjLB1 = 0;
            if ( ii == ii_Start ) {
                jjLB1 = jj_Start;
            }
            for ( int jj = jjLB1; jj < 32; jj++) {
                {
                int kkLB2 = 0;
                if ( ii == ii_Start && jj == jj_Start ) {
                    kkLB2 = kk_Start;
                }
                for ( int kk = kkLB2; kk < 32; kk++) {
                    {
                    int iLB3 = ii;
                    if ( ii == ii_Start && jj == jj_Start && kk == kk_Start ) {
                        iLB3 = i_Start;
                    }
                    for ( int i = iLB3; i < (ii + 8); i++) {
                        {
                        int jLB4 = jj;
                        if ( ii == ii_Start && jj == jj_Start && kk == kk_Start && i == i_Start ) {
                            jLB4 = j_Start;
                        }
                        for ( int j = jLB4; j < (jj + 8); j++) {
                            {
                            int kLB5 = kk;
                            if ( ii == ii_Start && jj == jj_Start && kk == kk_Start && i == i_Start && j == j_Start ) {
                                kLB5 = k_Start;
                            }
                            for ( int k = kLB5; k < (kk + 8); k++) {
                                if (cntStart == true) cnt++;
                                if (cntStart == true) {
                                    cnt++;
                                    if ( calAddrpB_addr1( ii, jj, kk, i, j, k) == calAddrpB_addr1(ii_Start, jj_Start, kk_Start, i_Start, j_Start, k_Start)) {
                                        rtHistoCal(cnt, 1);
                                        goto EndSample;
                                    }
                                }
                                cntStart = true;
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
    ref_pC_addr2();
    ref_pA_addr0();
    ref_pB_addr1();
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
 /* Analyze function: mm_tile */ 
