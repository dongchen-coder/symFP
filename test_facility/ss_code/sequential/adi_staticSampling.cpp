
 /* Start to analysis array index
Array index info: Total number of references: 34
p.addr (((i * 1024) + j) - 1)
u.addr ((j * 1024) + i)
u.addr (((j * 1024) + i) + 1)
v.addr (0 + i)
p.addr ((i * 1024) + 0)
v.addr (0 + i)
q.addr ((i * 1024) + 0)
v.addr (1047552 + i)
p.addr ((i * 1024) + j)
v.addr (((j + 1) * 1024) + i)
q.addr ((i * 1024) + j)
p.addr ((i * 1024) + j)
u.addr ((i * 1024) + 0)
u.addr (((j * 1024) + i) - 1)
q.addr ((i * 1024) + 0)
p.addr (((i * 1024) + j) - 1)
p.addr ((i * 1024) + j)
v.addr (((i - 1) * 1024) + j)
q.addr (((i * 1024) + j) - 1)
p.addr (((i * 1024) + j) - 1)
q.addr ((i * 1024) + j)
u.addr (((i * 1024) + j) + 1)
q.addr ((i * 1024) + j)
u.addr ((i * 1024) + j)
v.addr ((j * 1024) + i)
p.addr ((i * 1024) + 0)
u.addr ((i * 1024) + 0)
v.addr ((i * 1024) + j)
v.addr (((i + 1) * 1024) + j)
q.addr (((i * 1024) + j) - 1)
p.addr (((i * 1024) + j) - 1)
q.addr ((i * 1024) + j)
u.addr (((i * 1024) + 1024) - 1)
p.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %p
double* %q
double* %v
double* %u

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (1, 10)
--Loop inc: (t + 1)
--Loop predicate: <=
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------array access v.addr (0 + i)
------array access p.addr ((i * 1024) + 0)
------array access v.addr (0 + i)
------array access q.addr ((i * 1024) + 0)
------j
------Loop Bound: (1, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access p.addr (((i * 1024) + j) - 1)
--------array access p.addr ((i * 1024) + j)
--------array access u.addr (((j * 1024) + i) - 1)
--------array access u.addr ((j * 1024) + i)
--------array access u.addr (((j * 1024) + i) + 1)
--------array access q.addr (((i * 1024) + j) - 1)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access q.addr ((i * 1024) + j)
------array access v.addr (1047552 + i)
------j
------Loop Bound: (1022, 1)
------Loop inc: (j + -1)
------Loop predicate: >=
--------array access p.addr ((i * 1024) + j)
--------array access v.addr (((j + 1) * 1024) + i)
--------array access q.addr ((i * 1024) + j)
--------array access v.addr ((j * 1024) + i)
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------array access u.addr ((i * 1024) + 0)
------array access p.addr ((i * 1024) + 0)
------array access u.addr ((i * 1024) + 0)
------array access q.addr ((i * 1024) + 0)
------j
------Loop Bound: (1, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access p.addr (((i * 1024) + j) - 1)
--------array access p.addr ((i * 1024) + j)
--------array access v.addr (((i - 1) * 1024) + j)
--------array access v.addr ((i * 1024) + j)
--------array access v.addr (((i + 1) * 1024) + j)
--------array access q.addr (((i * 1024) + j) - 1)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access q.addr ((i * 1024) + j)
------array access u.addr (((i * 1024) + 1024) - 1)
------j
------Loop Bound: (1022, 1)
------Loop inc: (j + -1)
------Loop predicate: >=
--------array access p.addr ((i * 1024) + j)
--------array access u.addr (((i * 1024) + j) + 1)
--------array access q.addr ((i * 1024) + j)
--------array access u.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 0
------Sample number: 25
--------Sample number: 1305
--------Sample number: 1305
------Sample number: 25
--------Sample number: 1305
--------Sample number: 1305
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
/* v_addr (0 + i) 0 */
int calAddrv_addr0( int t, int i) {
    int result = ((0 + i)) * 8 / 64;
    return result;
}
/* p_addr ((i * 1024) + 0) 1 */
int calAddrp_addr1( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* v_addr (0 + i) 2 */
int calAddrv_addr2( int t, int i) {
    int result = ((0 + i)) * 8 / 64;
    return result;
}
/* q_addr ((i * 1024) + 0) 3 */
int calAddrq_addr3( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* p_addr (((i * 1024) + j) - 1) 4 */
int calAddrp_addr4( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* p_addr ((i * 1024) + j) 5 */
int calAddrp_addr5( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* u_addr (((j * 1024) + i) - 1) 6 */
int calAddru_addr6( int t, int i, int j) {
    int result = ((((j * 1024) + i) - 1)) * 8 / 64;
    return result;
}
/* u_addr ((j * 1024) + i) 7 */
int calAddru_addr7( int t, int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 64;
    return result;
}
/* u_addr (((j * 1024) + i) + 1) 8 */
int calAddru_addr8( int t, int i, int j) {
    int result = ((((j * 1024) + i) + 1)) * 8 / 64;
    return result;
}
/* q_addr (((i * 1024) + j) - 1) 9 */
int calAddrq_addr9( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* p_addr (((i * 1024) + j) - 1) 10 */
int calAddrp_addr10( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* q_addr ((i * 1024) + j) 11 */
int calAddrq_addr11( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* v_addr (1047552 + i) 12 */
int calAddrv_addr12( int t, int i) {
    int result = ((1047552 + i)) * 8 / 64;
    return result;
}
/* p_addr ((i * 1024) + j) 13 */
int calAddrp_addr13( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* v_addr (((j + 1) * 1024) + i) 14 */
int calAddrv_addr14( int t, int i, int j) {
    int result = ((((j + 1) * 1024) + i)) * 8 / 64;
    return result;
}
/* q_addr ((i * 1024) + j) 15 */
int calAddrq_addr15( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* v_addr ((j * 1024) + i) 16 */
int calAddrv_addr16( int t, int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 64;
    return result;
}
/* u_addr ((i * 1024) + 0) 17 */
int calAddru_addr17( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* p_addr ((i * 1024) + 0) 18 */
int calAddrp_addr18( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* u_addr ((i * 1024) + 0) 19 */
int calAddru_addr19( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* q_addr ((i * 1024) + 0) 20 */
int calAddrq_addr20( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
/* p_addr (((i * 1024) + j) - 1) 21 */
int calAddrp_addr21( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* p_addr ((i * 1024) + j) 22 */
int calAddrp_addr22( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* v_addr (((i - 1) * 1024) + j) 23 */
int calAddrv_addr23( int t, int i, int j) {
    int result = ((((i - 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* v_addr ((i * 1024) + j) 24 */
int calAddrv_addr24( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* v_addr (((i + 1) * 1024) + j) 25 */
int calAddrv_addr25( int t, int i, int j) {
    int result = ((((i + 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* q_addr (((i * 1024) + j) - 1) 26 */
int calAddrq_addr26( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* p_addr (((i * 1024) + j) - 1) 27 */
int calAddrp_addr27( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* q_addr ((i * 1024) + j) 28 */
int calAddrq_addr28( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* u_addr (((i * 1024) + 1024) - 1) 29 */
int calAddru_addr29( int t, int i) {
    int result = ((((i * 1024) + 1024) - 1)) * 8 / 64;
    return result;
}
/* p_addr ((i * 1024) + j) 30 */
int calAddrp_addr30( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* u_addr (((i * 1024) + j) + 1) 31 */
int calAddru_addr31( int t, int i, int j) {
    int result = ((((i * 1024) + j) + 1)) * 8 / 64;
    return result;
}
/* q_addr ((i * 1024) + j) 32 */
int calAddrq_addr32( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* u_addr ((i * 1024) + j) 33 */
int calAddru_addr33( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_p_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr4 = -1;
    uint64_t prev_t_Start_p_addr4 = -1;
    uint64_t prev_t_End_p_addr4 = -1;
    uint64_t prev_i_Start_p_addr4 = -1;
    uint64_t prev_i_End_p_addr4 = -1;
    uint64_t prev_j_Start_p_addr4 = -1;
    uint64_t prev_j_End_p_addr4 = -1;
    uint64_t prev_cnt_p_addr5 = -1;
    uint64_t prev_t_Start_p_addr5 = -1;
    uint64_t prev_t_End_p_addr5 = -1;
    uint64_t prev_i_Start_p_addr5 = -1;
    uint64_t prev_i_End_p_addr5 = -1;
    uint64_t prev_j_Start_p_addr5 = -1;
    uint64_t prev_j_End_p_addr5 = -1;
    uint64_t prev_cnt_p_addr10 = -1;
    uint64_t prev_t_Start_p_addr10 = -1;
    uint64_t prev_t_End_p_addr10 = -1;
    uint64_t prev_i_Start_p_addr10 = -1;
    uint64_t prev_i_End_p_addr10 = -1;
    uint64_t prev_j_Start_p_addr10 = -1;
    uint64_t prev_j_End_p_addr10 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr4 != -1) {
            if ( calAddrp_addr4( t_Start - prev_t_Start_p_addr4 + prev_t_End_p_addr4, i_Start - prev_i_Start_p_addr4 + prev_i_End_p_addr4, j_Start - prev_j_Start_p_addr4 + prev_j_End_p_addr4) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr4, 4);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr5 != -1) {
            if ( calAddrp_addr5( t_Start - prev_t_Start_p_addr5 + prev_t_End_p_addr5, i_Start - prev_i_Start_p_addr5 + prev_i_End_p_addr5, j_Start - prev_j_Start_p_addr5 + prev_j_End_p_addr5) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr5, 5);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr10 != -1) {
            if ( calAddrp_addr10( t_Start - prev_t_Start_p_addr10 + prev_t_End_p_addr10, i_Start - prev_i_Start_p_addr10 + prev_i_End_p_addr10, j_Start - prev_j_Start_p_addr10 + prev_j_End_p_addr10) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr10, 10);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 4);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 4);
                            prev_cnt_p_addr4 = cnt;
                            prev_t_Start_p_addr4 = t_Start;
                            prev_t_End_p_addr4 = t;
                            prev_i_Start_p_addr4 = i_Start;
                            prev_i_End_p_addr4 = i;
                            prev_j_Start_p_addr4 = j_Start;
                            prev_j_End_p_addr4 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 4);
                            prev_cnt_p_addr5 = cnt;
                            prev_t_Start_p_addr5 = t_Start;
                            prev_t_End_p_addr5 = t;
                            prev_i_Start_p_addr5 = i_Start;
                            prev_i_End_p_addr5 = i;
                            prev_j_Start_p_addr5 = j_Start;
                            prev_j_End_p_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 4);
                            prev_cnt_p_addr10 = cnt;
                            prev_t_Start_p_addr10 = t_Start;
                            prev_t_End_p_addr10 = t;
                            prev_i_Start_p_addr10 = i_Start;
                            prev_i_End_p_addr10 = i;
                            prev_j_Start_p_addr10 = j_Start;
                            prev_j_End_p_addr10 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
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
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 4);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 4);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 4);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 4);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
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
void ref_u_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr6 = -1;
    uint64_t prev_t_Start_u_addr6 = -1;
    uint64_t prev_t_End_u_addr6 = -1;
    uint64_t prev_i_Start_u_addr6 = -1;
    uint64_t prev_i_End_u_addr6 = -1;
    uint64_t prev_j_Start_u_addr6 = -1;
    uint64_t prev_j_End_u_addr6 = -1;
    uint64_t prev_cnt_u_addr7 = -1;
    uint64_t prev_t_Start_u_addr7 = -1;
    uint64_t prev_t_End_u_addr7 = -1;
    uint64_t prev_i_Start_u_addr7 = -1;
    uint64_t prev_i_End_u_addr7 = -1;
    uint64_t prev_j_Start_u_addr7 = -1;
    uint64_t prev_j_End_u_addr7 = -1;
    uint64_t prev_cnt_u_addr8 = -1;
    uint64_t prev_t_Start_u_addr8 = -1;
    uint64_t prev_t_End_u_addr8 = -1;
    uint64_t prev_i_Start_u_addr8 = -1;
    uint64_t prev_i_End_u_addr8 = -1;
    uint64_t prev_j_Start_u_addr8 = -1;
    uint64_t prev_j_End_u_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr6 != -1) {
            if ( calAddru_addr6( t_Start - prev_t_Start_u_addr6 + prev_t_End_u_addr6, i_Start - prev_i_Start_u_addr6 + prev_i_End_u_addr6, j_Start - prev_j_Start_u_addr6 + prev_j_End_u_addr6) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr6, 6);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr7 != -1) {
            if ( calAddru_addr7( t_Start - prev_t_Start_u_addr7 + prev_t_End_u_addr7, i_Start - prev_i_Start_u_addr7 + prev_i_End_u_addr7, j_Start - prev_j_Start_u_addr7 + prev_j_End_u_addr7) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr7, 7);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr8 != -1) {
            if ( calAddru_addr8( t_Start - prev_t_Start_u_addr8 + prev_t_End_u_addr8, i_Start - prev_i_Start_u_addr8 + prev_i_End_u_addr8, j_Start - prev_j_Start_u_addr8 + prev_j_End_u_addr8) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr8, 8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 7);
                            prev_cnt_u_addr6 = cnt;
                            prev_t_Start_u_addr6 = t_Start;
                            prev_t_End_u_addr6 = t;
                            prev_i_Start_u_addr6 = i_Start;
                            prev_i_End_u_addr6 = i;
                            prev_j_Start_u_addr6 = j_Start;
                            prev_j_End_u_addr6 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 7);
                            prev_cnt_u_addr7 = cnt;
                            prev_t_Start_u_addr7 = t_Start;
                            prev_t_End_u_addr7 = t;
                            prev_i_Start_u_addr7 = i_Start;
                            prev_i_End_u_addr7 = i;
                            prev_j_Start_u_addr7 = j_Start;
                            prev_j_End_u_addr7 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 7);
                            prev_cnt_u_addr8 = cnt;
                            prev_t_Start_u_addr8 = t_Start;
                            prev_t_End_u_addr8 = t;
                            prev_i_Start_u_addr8 = i_Start;
                            prev_i_End_u_addr8 = i;
                            prev_j_Start_u_addr8 = j_Start;
                            prev_j_End_u_addr8 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 7);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 7);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr29( t, i) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 7);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 7);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 7);
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
void ref_u_addr8() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr6 = -1;
    uint64_t prev_t_Start_u_addr6 = -1;
    uint64_t prev_t_End_u_addr6 = -1;
    uint64_t prev_i_Start_u_addr6 = -1;
    uint64_t prev_i_End_u_addr6 = -1;
    uint64_t prev_j_Start_u_addr6 = -1;
    uint64_t prev_j_End_u_addr6 = -1;
    uint64_t prev_cnt_u_addr7 = -1;
    uint64_t prev_t_Start_u_addr7 = -1;
    uint64_t prev_t_End_u_addr7 = -1;
    uint64_t prev_i_Start_u_addr7 = -1;
    uint64_t prev_i_End_u_addr7 = -1;
    uint64_t prev_j_Start_u_addr7 = -1;
    uint64_t prev_j_End_u_addr7 = -1;
    uint64_t prev_cnt_u_addr8 = -1;
    uint64_t prev_t_Start_u_addr8 = -1;
    uint64_t prev_t_End_u_addr8 = -1;
    uint64_t prev_i_Start_u_addr8 = -1;
    uint64_t prev_i_End_u_addr8 = -1;
    uint64_t prev_j_Start_u_addr8 = -1;
    uint64_t prev_j_End_u_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr6 != -1) {
            if ( calAddru_addr6( t_Start - prev_t_Start_u_addr6 + prev_t_End_u_addr6, i_Start - prev_i_Start_u_addr6 + prev_i_End_u_addr6, j_Start - prev_j_Start_u_addr6 + prev_j_End_u_addr6) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr6, 6);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr7 != -1) {
            if ( calAddru_addr7( t_Start - prev_t_Start_u_addr7 + prev_t_End_u_addr7, i_Start - prev_i_Start_u_addr7 + prev_i_End_u_addr7, j_Start - prev_j_Start_u_addr7 + prev_j_End_u_addr7) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr7, 7);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr8 != -1) {
            if ( calAddru_addr8( t_Start - prev_t_Start_u_addr8 + prev_t_End_u_addr8, i_Start - prev_i_Start_u_addr8 + prev_i_End_u_addr8, j_Start - prev_j_Start_u_addr8 + prev_j_End_u_addr8) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr8, 8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 8);
                            prev_cnt_u_addr6 = cnt;
                            prev_t_Start_u_addr6 = t_Start;
                            prev_t_End_u_addr6 = t;
                            prev_i_Start_u_addr6 = i_Start;
                            prev_i_End_u_addr6 = i;
                            prev_j_Start_u_addr6 = j_Start;
                            prev_j_End_u_addr6 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 8);
                            prev_cnt_u_addr7 = cnt;
                            prev_t_Start_u_addr7 = t_Start;
                            prev_t_End_u_addr7 = t;
                            prev_i_Start_u_addr7 = i_Start;
                            prev_i_End_u_addr7 = i;
                            prev_j_Start_u_addr7 = j_Start;
                            prev_j_End_u_addr7 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 8);
                            prev_cnt_u_addr8 = cnt;
                            prev_t_Start_u_addr8 = t_Start;
                            prev_t_End_u_addr8 = t;
                            prev_i_Start_u_addr8 = i_Start;
                            prev_i_End_u_addr8 = i;
                            prev_j_Start_u_addr8 = j_Start;
                            prev_j_End_u_addr8 = j;
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
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 8);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 8);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr29( t, i) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 8);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 8);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 8);
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
void ref_v_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr0 = -1;
    uint64_t prev_t_Start_v_addr0 = -1;
    uint64_t prev_t_End_v_addr0 = -1;
    uint64_t prev_i_Start_v_addr0 = -1;
    uint64_t prev_i_End_v_addr0 = -1;
    uint64_t prev_cnt_v_addr2 = -1;
    uint64_t prev_t_Start_v_addr2 = -1;
    uint64_t prev_t_End_v_addr2 = -1;
    uint64_t prev_i_Start_v_addr2 = -1;
    uint64_t prev_i_End_v_addr2 = -1;
    uint64_t prev_cnt_v_addr12 = -1;
    uint64_t prev_t_Start_v_addr12 = -1;
    uint64_t prev_t_End_v_addr12 = -1;
    uint64_t prev_i_Start_v_addr12 = -1;
    uint64_t prev_i_End_v_addr12 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr0 != -1) {
            if ( calAddrv_addr0( t_Start - prev_t_Start_v_addr0 + prev_t_End_v_addr0, i_Start - prev_i_Start_v_addr0 + prev_i_End_v_addr0) == calAddrv_addr0(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr0, 0);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr2 != -1) {
            if ( calAddrv_addr2( t_Start - prev_t_Start_v_addr2 + prev_t_End_v_addr2, i_Start - prev_i_Start_v_addr2 + prev_i_End_v_addr2) == calAddrv_addr0(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr2, 2);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr12 != -1) {
            if ( calAddrv_addr12( t_Start - prev_t_Start_v_addr12 + prev_t_End_v_addr12, i_Start - prev_i_Start_v_addr12 + prev_i_End_v_addr12) == calAddrv_addr0(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr12, 12);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt, 0);
                        prev_cnt_v_addr0 = cnt;
                        prev_t_Start_v_addr0 = t_Start;
                        prev_t_End_v_addr0 = t;
                        prev_i_Start_v_addr0 = i_Start;
                        prev_i_End_v_addr0 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt, 0);
                        prev_cnt_v_addr2 = cnt;
                        prev_t_Start_v_addr2 = t_Start;
                        prev_t_End_v_addr2 = t;
                        prev_i_Start_v_addr2 = i_Start;
                        prev_i_End_v_addr2 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr12( t, i) == calAddrv_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt, 0);
                        prev_cnt_v_addr12 = cnt;
                        prev_t_Start_v_addr12 = t_Start;
                        prev_t_End_v_addr12 = t;
                        prev_i_Start_v_addr12 = i_Start;
                        prev_i_End_v_addr12 = i;
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt, 0);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt, 0);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt, 0);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt, 0);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt, 0);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
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
void ref_p_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr1 = -1;
    uint64_t prev_t_Start_p_addr1 = -1;
    uint64_t prev_t_End_p_addr1 = -1;
    uint64_t prev_i_Start_p_addr1 = -1;
    uint64_t prev_i_End_p_addr1 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr1 != -1) {
            if ( calAddrp_addr1( t_Start - prev_t_Start_p_addr1 + prev_t_End_p_addr1, i_Start - prev_i_Start_p_addr1 + prev_i_End_p_addr1) == calAddrp_addr1(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_p_addr1, 1);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt, 1);
                        prev_cnt_p_addr1 = cnt;
                        prev_t_Start_p_addr1 = t_Start;
                        prev_t_End_p_addr1 = t;
                        prev_i_Start_p_addr1 = i_Start;
                        prev_i_End_p_addr1 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt, 1);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt, 1);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt, 1);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt, 1);
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt, 1);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt, 1);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt, 1);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt, 1);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt, 1);
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
void ref_v_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr0 = -1;
    uint64_t prev_t_Start_v_addr0 = -1;
    uint64_t prev_t_End_v_addr0 = -1;
    uint64_t prev_i_Start_v_addr0 = -1;
    uint64_t prev_i_End_v_addr0 = -1;
    uint64_t prev_cnt_v_addr2 = -1;
    uint64_t prev_t_Start_v_addr2 = -1;
    uint64_t prev_t_End_v_addr2 = -1;
    uint64_t prev_i_Start_v_addr2 = -1;
    uint64_t prev_i_End_v_addr2 = -1;
    uint64_t prev_cnt_v_addr12 = -1;
    uint64_t prev_t_Start_v_addr12 = -1;
    uint64_t prev_t_End_v_addr12 = -1;
    uint64_t prev_i_Start_v_addr12 = -1;
    uint64_t prev_i_End_v_addr12 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr0 != -1) {
            if ( calAddrv_addr0( t_Start - prev_t_Start_v_addr0 + prev_t_End_v_addr0, i_Start - prev_i_Start_v_addr0 + prev_i_End_v_addr0) == calAddrv_addr2(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr0, 0);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr2 != -1) {
            if ( calAddrv_addr2( t_Start - prev_t_Start_v_addr2 + prev_t_End_v_addr2, i_Start - prev_i_Start_v_addr2 + prev_i_End_v_addr2) == calAddrv_addr2(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr2, 2);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr12 != -1) {
            if ( calAddrv_addr12( t_Start - prev_t_Start_v_addr12 + prev_t_End_v_addr12, i_Start - prev_i_Start_v_addr12 + prev_i_End_v_addr12) == calAddrv_addr2(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr12, 12);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt, 2);
                        prev_cnt_v_addr0 = cnt;
                        prev_t_Start_v_addr0 = t_Start;
                        prev_t_End_v_addr0 = t;
                        prev_i_Start_v_addr0 = i_Start;
                        prev_i_End_v_addr0 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt, 2);
                        prev_cnt_v_addr2 = cnt;
                        prev_t_Start_v_addr2 = t_Start;
                        prev_t_End_v_addr2 = t;
                        prev_i_Start_v_addr2 = i_Start;
                        prev_i_End_v_addr2 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr12( t, i) == calAddrv_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt, 2);
                        prev_cnt_v_addr12 = cnt;
                        prev_t_Start_v_addr12 = t_Start;
                        prev_t_End_v_addr12 = t;
                        prev_i_Start_v_addr12 = i_Start;
                        prev_i_End_v_addr12 = i;
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt, 2);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt, 2);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt, 2);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt, 2);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt, 2);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
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
void ref_q_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr3 = -1;
    uint64_t prev_t_Start_q_addr3 = -1;
    uint64_t prev_t_End_q_addr3 = -1;
    uint64_t prev_i_Start_q_addr3 = -1;
    uint64_t prev_i_End_q_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr3 != -1) {
            if ( calAddrq_addr3( t_Start - prev_t_Start_q_addr3 + prev_t_End_q_addr3, i_Start - prev_i_Start_q_addr3 + prev_i_End_q_addr3) == calAddrq_addr3(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_q_addr3, 3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt, 3);
                        prev_cnt_q_addr3 = cnt;
                        prev_t_Start_q_addr3 = t_Start;
                        prev_t_End_q_addr3 = t;
                        prev_i_Start_q_addr3 = i_Start;
                        prev_i_End_q_addr3 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt, 3);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt, 3);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt, 3);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt, 3);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt, 3);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt, 3);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt, 3);
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
void ref_v_addr12() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr0 = -1;
    uint64_t prev_t_Start_v_addr0 = -1;
    uint64_t prev_t_End_v_addr0 = -1;
    uint64_t prev_i_Start_v_addr0 = -1;
    uint64_t prev_i_End_v_addr0 = -1;
    uint64_t prev_cnt_v_addr2 = -1;
    uint64_t prev_t_Start_v_addr2 = -1;
    uint64_t prev_t_End_v_addr2 = -1;
    uint64_t prev_i_Start_v_addr2 = -1;
    uint64_t prev_i_End_v_addr2 = -1;
    uint64_t prev_cnt_v_addr12 = -1;
    uint64_t prev_t_Start_v_addr12 = -1;
    uint64_t prev_t_End_v_addr12 = -1;
    uint64_t prev_i_Start_v_addr12 = -1;
    uint64_t prev_i_End_v_addr12 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr0 != -1) {
            if ( calAddrv_addr0( t_Start - prev_t_Start_v_addr0 + prev_t_End_v_addr0, i_Start - prev_i_Start_v_addr0 + prev_i_End_v_addr0) == calAddrv_addr12(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr0, 0);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr2 != -1) {
            if ( calAddrv_addr2( t_Start - prev_t_Start_v_addr2 + prev_t_End_v_addr2, i_Start - prev_i_Start_v_addr2 + prev_i_End_v_addr2) == calAddrv_addr12(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr2, 2);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr12 != -1) {
            if ( calAddrv_addr12( t_Start - prev_t_Start_v_addr12 + prev_t_End_v_addr12, i_Start - prev_i_Start_v_addr12 + prev_i_End_v_addr12) == calAddrv_addr12(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr12, 12);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr12(t_Start, i_Start)) {
                        rtHistoCal(cnt, 12);
                        prev_cnt_v_addr0 = cnt;
                        prev_t_Start_v_addr0 = t_Start;
                        prev_t_End_v_addr0 = t;
                        prev_i_Start_v_addr0 = i_Start;
                        prev_i_End_v_addr0 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr12(t_Start, i_Start)) {
                        rtHistoCal(cnt, 12);
                        prev_cnt_v_addr2 = cnt;
                        prev_t_Start_v_addr2 = t_Start;
                        prev_t_End_v_addr2 = t;
                        prev_i_Start_v_addr2 = i_Start;
                        prev_i_End_v_addr2 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr12( t, i) == calAddrv_addr12(t_Start, i_Start)) {
                        rtHistoCal(cnt, 12);
                        prev_cnt_v_addr12 = cnt;
                        prev_t_Start_v_addr12 = t_Start;
                        prev_t_End_v_addr12 = t;
                        prev_i_Start_v_addr12 = i_Start;
                        prev_i_End_v_addr12 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                            rtHistoCal(cnt, 12);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                            rtHistoCal(cnt, 12);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                            rtHistoCal(cnt, 12);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                            rtHistoCal(cnt, 12);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr12(t_Start, i_Start)) {
                            rtHistoCal(cnt, 12);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
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
void ref_p_addr13() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr13 = -1;
    uint64_t prev_t_Start_p_addr13 = -1;
    uint64_t prev_t_End_p_addr13 = -1;
    uint64_t prev_i_Start_p_addr13 = -1;
    uint64_t prev_i_End_p_addr13 = -1;
    uint64_t prev_j_Start_p_addr13 = -1;
    uint64_t prev_j_End_p_addr13 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr13 != -1) {
            if ( calAddrp_addr13( t_Start - prev_t_Start_p_addr13 + prev_t_End_p_addr13, i_Start - prev_i_Start_p_addr13 + prev_i_End_p_addr13, j_Start - prev_j_Start_p_addr13 + prev_j_End_p_addr13) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr13, 13);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 13);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 13);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 13);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 13);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 13);
                            prev_cnt_p_addr13 = cnt;
                            prev_t_Start_p_addr13 = t_Start;
                            prev_t_End_p_addr13 = t;
                            prev_i_Start_p_addr13 = i_Start;
                            prev_i_End_p_addr13 = i;
                            prev_j_Start_p_addr13 = j_Start;
                            prev_j_End_p_addr13 = j;
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 13);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 13);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 13);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 13);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr13(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 13);
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
void ref_v_addr14() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr14 = -1;
    uint64_t prev_t_Start_v_addr14 = -1;
    uint64_t prev_t_End_v_addr14 = -1;
    uint64_t prev_i_Start_v_addr14 = -1;
    uint64_t prev_i_End_v_addr14 = -1;
    uint64_t prev_j_Start_v_addr14 = -1;
    uint64_t prev_j_End_v_addr14 = -1;
    uint64_t prev_cnt_v_addr16 = -1;
    uint64_t prev_t_Start_v_addr16 = -1;
    uint64_t prev_t_End_v_addr16 = -1;
    uint64_t prev_i_Start_v_addr16 = -1;
    uint64_t prev_i_End_v_addr16 = -1;
    uint64_t prev_j_Start_v_addr16 = -1;
    uint64_t prev_j_End_v_addr16 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr14 != -1) {
            if ( calAddrv_addr14( t_Start - prev_t_Start_v_addr14 + prev_t_End_v_addr14, i_Start - prev_i_Start_v_addr14 + prev_i_End_v_addr14, j_Start - prev_j_Start_v_addr14 + prev_j_End_v_addr14) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr14, 14);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr16 != -1) {
            if ( calAddrv_addr16( t_Start - prev_t_Start_v_addr16 + prev_t_End_v_addr16, i_Start - prev_i_Start_v_addr16 + prev_i_End_v_addr16, j_Start - prev_j_Start_v_addr16 + prev_j_End_v_addr16) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr16, 16);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 14);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 14);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr12( t, i) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 14);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 14);
                            prev_cnt_v_addr14 = cnt;
                            prev_t_Start_v_addr14 = t_Start;
                            prev_t_End_v_addr14 = t;
                            prev_i_Start_v_addr14 = i_Start;
                            prev_i_End_v_addr14 = i;
                            prev_j_Start_v_addr14 = j_Start;
                            prev_j_End_v_addr14 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 14);
                            prev_cnt_v_addr16 = cnt;
                            prev_t_Start_v_addr16 = t_Start;
                            prev_t_End_v_addr16 = t;
                            prev_i_Start_v_addr16 = i_Start;
                            prev_i_End_v_addr16 = i;
                            prev_j_Start_v_addr16 = j_Start;
                            prev_j_End_v_addr16 = j;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 14);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 14);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr14(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 14);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
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
void ref_q_addr15() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr15 = -1;
    uint64_t prev_t_Start_q_addr15 = -1;
    uint64_t prev_t_End_q_addr15 = -1;
    uint64_t prev_i_Start_q_addr15 = -1;
    uint64_t prev_i_End_q_addr15 = -1;
    uint64_t prev_j_Start_q_addr15 = -1;
    uint64_t prev_j_End_q_addr15 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr15 != -1) {
            if ( calAddrq_addr15( t_Start - prev_t_Start_q_addr15 + prev_t_End_q_addr15, i_Start - prev_i_Start_q_addr15 + prev_i_End_q_addr15, j_Start - prev_j_Start_q_addr15 + prev_j_End_q_addr15) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr15, 15);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 15);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 15);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 15);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 15);
                            prev_cnt_q_addr15 = cnt;
                            prev_t_Start_q_addr15 = t_Start;
                            prev_t_End_q_addr15 = t;
                            prev_i_Start_q_addr15 = i_Start;
                            prev_i_End_q_addr15 = i;
                            prev_j_Start_q_addr15 = j_Start;
                            prev_j_End_q_addr15 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 15);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 15);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 15);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr15(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 15);
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
void ref_p_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr4 = -1;
    uint64_t prev_t_Start_p_addr4 = -1;
    uint64_t prev_t_End_p_addr4 = -1;
    uint64_t prev_i_Start_p_addr4 = -1;
    uint64_t prev_i_End_p_addr4 = -1;
    uint64_t prev_j_Start_p_addr4 = -1;
    uint64_t prev_j_End_p_addr4 = -1;
    uint64_t prev_cnt_p_addr5 = -1;
    uint64_t prev_t_Start_p_addr5 = -1;
    uint64_t prev_t_End_p_addr5 = -1;
    uint64_t prev_i_Start_p_addr5 = -1;
    uint64_t prev_i_End_p_addr5 = -1;
    uint64_t prev_j_Start_p_addr5 = -1;
    uint64_t prev_j_End_p_addr5 = -1;
    uint64_t prev_cnt_p_addr10 = -1;
    uint64_t prev_t_Start_p_addr10 = -1;
    uint64_t prev_t_End_p_addr10 = -1;
    uint64_t prev_i_Start_p_addr10 = -1;
    uint64_t prev_i_End_p_addr10 = -1;
    uint64_t prev_j_Start_p_addr10 = -1;
    uint64_t prev_j_End_p_addr10 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr4 != -1) {
            if ( calAddrp_addr4( t_Start - prev_t_Start_p_addr4 + prev_t_End_p_addr4, i_Start - prev_i_Start_p_addr4 + prev_i_End_p_addr4, j_Start - prev_j_Start_p_addr4 + prev_j_End_p_addr4) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr4, 4);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr5 != -1) {
            if ( calAddrp_addr5( t_Start - prev_t_Start_p_addr5 + prev_t_End_p_addr5, i_Start - prev_i_Start_p_addr5 + prev_i_End_p_addr5, j_Start - prev_j_Start_p_addr5 + prev_j_End_p_addr5) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr5, 5);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr10 != -1) {
            if ( calAddrp_addr10( t_Start - prev_t_Start_p_addr10 + prev_t_End_p_addr10, i_Start - prev_i_Start_p_addr10 + prev_i_End_p_addr10, j_Start - prev_j_Start_p_addr10 + prev_j_End_p_addr10) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr10, 10);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 5);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
                            prev_cnt_p_addr4 = cnt;
                            prev_t_Start_p_addr4 = t_Start;
                            prev_t_End_p_addr4 = t;
                            prev_i_Start_p_addr4 = i_Start;
                            prev_i_End_p_addr4 = i;
                            prev_j_Start_p_addr4 = j_Start;
                            prev_j_End_p_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
                            prev_cnt_p_addr5 = cnt;
                            prev_t_Start_p_addr5 = t_Start;
                            prev_t_End_p_addr5 = t;
                            prev_i_Start_p_addr5 = i_Start;
                            prev_i_End_p_addr5 = i;
                            prev_j_Start_p_addr5 = j_Start;
                            prev_j_End_p_addr5 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
                            prev_cnt_p_addr10 = cnt;
                            prev_t_Start_p_addr10 = t_Start;
                            prev_t_End_p_addr10 = t;
                            prev_i_Start_p_addr10 = i_Start;
                            prev_i_End_p_addr10 = i;
                            prev_j_Start_p_addr10 = j_Start;
                            prev_j_End_p_addr10 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 5);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 5);
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
void ref_u_addr17() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr17 = -1;
    uint64_t prev_t_Start_u_addr17 = -1;
    uint64_t prev_t_End_u_addr17 = -1;
    uint64_t prev_i_Start_u_addr17 = -1;
    uint64_t prev_i_End_u_addr17 = -1;
    uint64_t prev_cnt_u_addr19 = -1;
    uint64_t prev_t_Start_u_addr19 = -1;
    uint64_t prev_t_End_u_addr19 = -1;
    uint64_t prev_i_Start_u_addr19 = -1;
    uint64_t prev_i_End_u_addr19 = -1;
    uint64_t prev_cnt_u_addr29 = -1;
    uint64_t prev_t_Start_u_addr29 = -1;
    uint64_t prev_t_End_u_addr29 = -1;
    uint64_t prev_i_Start_u_addr29 = -1;
    uint64_t prev_i_End_u_addr29 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr17 != -1) {
            if ( calAddru_addr17( t_Start - prev_t_Start_u_addr17 + prev_t_End_u_addr17, i_Start - prev_i_Start_u_addr17 + prev_i_End_u_addr17) == calAddru_addr17(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr17, 17);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr19 != -1) {
            if ( calAddru_addr19( t_Start - prev_t_Start_u_addr19 + prev_t_End_u_addr19, i_Start - prev_i_Start_u_addr19 + prev_i_End_u_addr19) == calAddru_addr17(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr19, 19);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr29 != -1) {
            if ( calAddru_addr29( t_Start - prev_t_Start_u_addr29 + prev_t_End_u_addr29, i_Start - prev_i_Start_u_addr29 + prev_i_End_u_addr29) == calAddru_addr17(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr29, 29);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                            rtHistoCal(cnt, 17);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                            rtHistoCal(cnt, 17);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                            rtHistoCal(cnt, 17);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr17(t_Start, i_Start)) {
                        rtHistoCal(cnt, 17);
                        prev_cnt_u_addr17 = cnt;
                        prev_t_Start_u_addr17 = t_Start;
                        prev_t_End_u_addr17 = t;
                        prev_i_Start_u_addr17 = i_Start;
                        prev_i_End_u_addr17 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr17(t_Start, i_Start)) {
                        rtHistoCal(cnt, 17);
                        prev_cnt_u_addr19 = cnt;
                        prev_t_Start_u_addr19 = t_Start;
                        prev_t_End_u_addr19 = t;
                        prev_i_Start_u_addr19 = i_Start;
                        prev_i_End_u_addr19 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr29( t, i) == calAddru_addr17(t_Start, i_Start)) {
                        rtHistoCal(cnt, 17);
                        prev_cnt_u_addr29 = cnt;
                        prev_t_Start_u_addr29 = t_Start;
                        prev_t_End_u_addr29 = t;
                        prev_i_Start_u_addr29 = i_Start;
                        prev_i_End_u_addr29 = i;
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                            rtHistoCal(cnt, 17);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr17(t_Start, i_Start)) {
                            rtHistoCal(cnt, 17);
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
void ref_u_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr6 = -1;
    uint64_t prev_t_Start_u_addr6 = -1;
    uint64_t prev_t_End_u_addr6 = -1;
    uint64_t prev_i_Start_u_addr6 = -1;
    uint64_t prev_i_End_u_addr6 = -1;
    uint64_t prev_j_Start_u_addr6 = -1;
    uint64_t prev_j_End_u_addr6 = -1;
    uint64_t prev_cnt_u_addr7 = -1;
    uint64_t prev_t_Start_u_addr7 = -1;
    uint64_t prev_t_End_u_addr7 = -1;
    uint64_t prev_i_Start_u_addr7 = -1;
    uint64_t prev_i_End_u_addr7 = -1;
    uint64_t prev_j_Start_u_addr7 = -1;
    uint64_t prev_j_End_u_addr7 = -1;
    uint64_t prev_cnt_u_addr8 = -1;
    uint64_t prev_t_Start_u_addr8 = -1;
    uint64_t prev_t_End_u_addr8 = -1;
    uint64_t prev_i_Start_u_addr8 = -1;
    uint64_t prev_i_End_u_addr8 = -1;
    uint64_t prev_j_Start_u_addr8 = -1;
    uint64_t prev_j_End_u_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr6 != -1) {
            if ( calAddru_addr6( t_Start - prev_t_Start_u_addr6 + prev_t_End_u_addr6, i_Start - prev_i_Start_u_addr6 + prev_i_End_u_addr6, j_Start - prev_j_Start_u_addr6 + prev_j_End_u_addr6) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr6, 6);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr7 != -1) {
            if ( calAddru_addr7( t_Start - prev_t_Start_u_addr7 + prev_t_End_u_addr7, i_Start - prev_i_Start_u_addr7 + prev_i_End_u_addr7, j_Start - prev_j_Start_u_addr7 + prev_j_End_u_addr7) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr7, 7);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr8 != -1) {
            if ( calAddru_addr8( t_Start - prev_t_Start_u_addr8 + prev_t_End_u_addr8, i_Start - prev_i_Start_u_addr8 + prev_i_End_u_addr8, j_Start - prev_j_Start_u_addr8 + prev_j_End_u_addr8) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr8, 8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 6);
                            prev_cnt_u_addr6 = cnt;
                            prev_t_Start_u_addr6 = t_Start;
                            prev_t_End_u_addr6 = t;
                            prev_i_Start_u_addr6 = i_Start;
                            prev_i_End_u_addr6 = i;
                            prev_j_Start_u_addr6 = j_Start;
                            prev_j_End_u_addr6 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 6);
                            prev_cnt_u_addr7 = cnt;
                            prev_t_Start_u_addr7 = t_Start;
                            prev_t_End_u_addr7 = t;
                            prev_i_Start_u_addr7 = i_Start;
                            prev_i_End_u_addr7 = i;
                            prev_j_Start_u_addr7 = j_Start;
                            prev_j_End_u_addr7 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 6);
                            prev_cnt_u_addr8 = cnt;
                            prev_t_Start_u_addr8 = t_Start;
                            prev_t_End_u_addr8 = t;
                            prev_i_Start_u_addr8 = i_Start;
                            prev_i_End_u_addr8 = i;
                            prev_j_Start_u_addr8 = j_Start;
                            prev_j_End_u_addr8 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 6);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 6);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr29( t, i) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 6);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 6);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
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
void ref_q_addr20() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr20 = -1;
    uint64_t prev_t_Start_q_addr20 = -1;
    uint64_t prev_t_End_q_addr20 = -1;
    uint64_t prev_i_Start_q_addr20 = -1;
    uint64_t prev_i_End_q_addr20 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr20 != -1) {
            if ( calAddrq_addr20( t_Start - prev_t_Start_q_addr20 + prev_t_End_q_addr20, i_Start - prev_i_Start_q_addr20 + prev_i_End_q_addr20) == calAddrq_addr20(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_q_addr20, 20);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr20(t_Start, i_Start)) {
                        rtHistoCal(cnt, 20);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                            rtHistoCal(cnt, 20);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                            rtHistoCal(cnt, 20);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                            rtHistoCal(cnt, 20);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr20(t_Start, i_Start)) {
                        rtHistoCal(cnt, 20);
                        prev_cnt_q_addr20 = cnt;
                        prev_t_Start_q_addr20 = t_Start;
                        prev_t_End_q_addr20 = t;
                        prev_i_Start_q_addr20 = i_Start;
                        prev_i_End_q_addr20 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                            rtHistoCal(cnt, 20);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                            rtHistoCal(cnt, 20);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr20(t_Start, i_Start)) {
                            rtHistoCal(cnt, 20);
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
void ref_p_addr21() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr21 = -1;
    uint64_t prev_t_Start_p_addr21 = -1;
    uint64_t prev_t_End_p_addr21 = -1;
    uint64_t prev_i_Start_p_addr21 = -1;
    uint64_t prev_i_End_p_addr21 = -1;
    uint64_t prev_j_Start_p_addr21 = -1;
    uint64_t prev_j_End_p_addr21 = -1;
    uint64_t prev_cnt_p_addr22 = -1;
    uint64_t prev_t_Start_p_addr22 = -1;
    uint64_t prev_t_End_p_addr22 = -1;
    uint64_t prev_i_Start_p_addr22 = -1;
    uint64_t prev_i_End_p_addr22 = -1;
    uint64_t prev_j_Start_p_addr22 = -1;
    uint64_t prev_j_End_p_addr22 = -1;
    uint64_t prev_cnt_p_addr27 = -1;
    uint64_t prev_t_Start_p_addr27 = -1;
    uint64_t prev_t_End_p_addr27 = -1;
    uint64_t prev_i_Start_p_addr27 = -1;
    uint64_t prev_i_End_p_addr27 = -1;
    uint64_t prev_j_Start_p_addr27 = -1;
    uint64_t prev_j_End_p_addr27 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr21 != -1) {
            if ( calAddrp_addr21( t_Start - prev_t_Start_p_addr21 + prev_t_End_p_addr21, i_Start - prev_i_Start_p_addr21 + prev_i_End_p_addr21, j_Start - prev_j_Start_p_addr21 + prev_j_End_p_addr21) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr21, 21);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr22 != -1) {
            if ( calAddrp_addr22( t_Start - prev_t_Start_p_addr22 + prev_t_End_p_addr22, i_Start - prev_i_Start_p_addr22 + prev_i_End_p_addr22, j_Start - prev_j_Start_p_addr22 + prev_j_End_p_addr22) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr22, 22);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr27 != -1) {
            if ( calAddrp_addr27( t_Start - prev_t_Start_p_addr27 + prev_t_End_p_addr27, i_Start - prev_i_Start_p_addr27 + prev_i_End_p_addr27, j_Start - prev_j_Start_p_addr27 + prev_j_End_p_addr27) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr27, 27);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 21);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 21);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 21);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 21);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 21);
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 21);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 21);
                            prev_cnt_p_addr21 = cnt;
                            prev_t_Start_p_addr21 = t_Start;
                            prev_t_End_p_addr21 = t;
                            prev_i_Start_p_addr21 = i_Start;
                            prev_i_End_p_addr21 = i;
                            prev_j_Start_p_addr21 = j_Start;
                            prev_j_End_p_addr21 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 21);
                            prev_cnt_p_addr22 = cnt;
                            prev_t_Start_p_addr22 = t_Start;
                            prev_t_End_p_addr22 = t;
                            prev_i_Start_p_addr22 = i_Start;
                            prev_i_End_p_addr22 = i;
                            prev_j_Start_p_addr22 = j_Start;
                            prev_j_End_p_addr22 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 21);
                            prev_cnt_p_addr27 = cnt;
                            prev_t_Start_p_addr27 = t_Start;
                            prev_t_End_p_addr27 = t;
                            prev_i_Start_p_addr27 = i_Start;
                            prev_i_End_p_addr27 = i;
                            prev_j_Start_p_addr27 = j_Start;
                            prev_j_End_p_addr27 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr21(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 21);
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
void ref_p_addr22() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr21 = -1;
    uint64_t prev_t_Start_p_addr21 = -1;
    uint64_t prev_t_End_p_addr21 = -1;
    uint64_t prev_i_Start_p_addr21 = -1;
    uint64_t prev_i_End_p_addr21 = -1;
    uint64_t prev_j_Start_p_addr21 = -1;
    uint64_t prev_j_End_p_addr21 = -1;
    uint64_t prev_cnt_p_addr22 = -1;
    uint64_t prev_t_Start_p_addr22 = -1;
    uint64_t prev_t_End_p_addr22 = -1;
    uint64_t prev_i_Start_p_addr22 = -1;
    uint64_t prev_i_End_p_addr22 = -1;
    uint64_t prev_j_Start_p_addr22 = -1;
    uint64_t prev_j_End_p_addr22 = -1;
    uint64_t prev_cnt_p_addr27 = -1;
    uint64_t prev_t_Start_p_addr27 = -1;
    uint64_t prev_t_End_p_addr27 = -1;
    uint64_t prev_i_Start_p_addr27 = -1;
    uint64_t prev_i_End_p_addr27 = -1;
    uint64_t prev_j_Start_p_addr27 = -1;
    uint64_t prev_j_End_p_addr27 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr21 != -1) {
            if ( calAddrp_addr21( t_Start - prev_t_Start_p_addr21 + prev_t_End_p_addr21, i_Start - prev_i_Start_p_addr21 + prev_i_End_p_addr21, j_Start - prev_j_Start_p_addr21 + prev_j_End_p_addr21) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr21, 21);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr22 != -1) {
            if ( calAddrp_addr22( t_Start - prev_t_Start_p_addr22 + prev_t_End_p_addr22, i_Start - prev_i_Start_p_addr22 + prev_i_End_p_addr22, j_Start - prev_j_Start_p_addr22 + prev_j_End_p_addr22) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr22, 22);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr27 != -1) {
            if ( calAddrp_addr27( t_Start - prev_t_Start_p_addr27 + prev_t_End_p_addr27, i_Start - prev_i_Start_p_addr27 + prev_i_End_p_addr27, j_Start - prev_j_Start_p_addr27 + prev_j_End_p_addr27) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr27, 27);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 22);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 22);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 22);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 22);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 22);
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 22);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 22);
                            prev_cnt_p_addr21 = cnt;
                            prev_t_Start_p_addr21 = t_Start;
                            prev_t_End_p_addr21 = t;
                            prev_i_Start_p_addr21 = i_Start;
                            prev_i_End_p_addr21 = i;
                            prev_j_Start_p_addr21 = j_Start;
                            prev_j_End_p_addr21 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 22);
                            prev_cnt_p_addr22 = cnt;
                            prev_t_Start_p_addr22 = t_Start;
                            prev_t_End_p_addr22 = t;
                            prev_i_Start_p_addr22 = i_Start;
                            prev_i_End_p_addr22 = i;
                            prev_j_Start_p_addr22 = j_Start;
                            prev_j_End_p_addr22 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 22);
                            prev_cnt_p_addr27 = cnt;
                            prev_t_Start_p_addr27 = t_Start;
                            prev_t_End_p_addr27 = t;
                            prev_i_Start_p_addr27 = i_Start;
                            prev_i_End_p_addr27 = i;
                            prev_j_Start_p_addr27 = j_Start;
                            prev_j_End_p_addr27 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr22(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 22);
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
void ref_v_addr23() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr23 = -1;
    uint64_t prev_t_Start_v_addr23 = -1;
    uint64_t prev_t_End_v_addr23 = -1;
    uint64_t prev_i_Start_v_addr23 = -1;
    uint64_t prev_i_End_v_addr23 = -1;
    uint64_t prev_j_Start_v_addr23 = -1;
    uint64_t prev_j_End_v_addr23 = -1;
    uint64_t prev_cnt_v_addr24 = -1;
    uint64_t prev_t_Start_v_addr24 = -1;
    uint64_t prev_t_End_v_addr24 = -1;
    uint64_t prev_i_Start_v_addr24 = -1;
    uint64_t prev_i_End_v_addr24 = -1;
    uint64_t prev_j_Start_v_addr24 = -1;
    uint64_t prev_j_End_v_addr24 = -1;
    uint64_t prev_cnt_v_addr25 = -1;
    uint64_t prev_t_Start_v_addr25 = -1;
    uint64_t prev_t_End_v_addr25 = -1;
    uint64_t prev_i_Start_v_addr25 = -1;
    uint64_t prev_i_End_v_addr25 = -1;
    uint64_t prev_j_Start_v_addr25 = -1;
    uint64_t prev_j_End_v_addr25 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr23 != -1) {
            if ( calAddrv_addr23( t_Start - prev_t_Start_v_addr23 + prev_t_End_v_addr23, i_Start - prev_i_Start_v_addr23 + prev_i_End_v_addr23, j_Start - prev_j_Start_v_addr23 + prev_j_End_v_addr23) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr23, 23);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr24 != -1) {
            if ( calAddrv_addr24( t_Start - prev_t_Start_v_addr24 + prev_t_End_v_addr24, i_Start - prev_i_Start_v_addr24 + prev_i_End_v_addr24, j_Start - prev_j_Start_v_addr24 + prev_j_End_v_addr24) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr24, 24);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr25 != -1) {
            if ( calAddrv_addr25( t_Start - prev_t_Start_v_addr25 + prev_t_End_v_addr25, i_Start - prev_i_Start_v_addr25 + prev_i_End_v_addr25, j_Start - prev_j_Start_v_addr25 + prev_j_End_v_addr25) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr25, 25);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 23);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 23);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr12( t, i) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 23);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 23);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 23);
                            goto EndSample;
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
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 23);
                            prev_cnt_v_addr23 = cnt;
                            prev_t_Start_v_addr23 = t_Start;
                            prev_t_End_v_addr23 = t;
                            prev_i_Start_v_addr23 = i_Start;
                            prev_i_End_v_addr23 = i;
                            prev_j_Start_v_addr23 = j_Start;
                            prev_j_End_v_addr23 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 23);
                            prev_cnt_v_addr24 = cnt;
                            prev_t_Start_v_addr24 = t_Start;
                            prev_t_End_v_addr24 = t;
                            prev_i_Start_v_addr24 = i_Start;
                            prev_i_End_v_addr24 = i;
                            prev_j_Start_v_addr24 = j_Start;
                            prev_j_End_v_addr24 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr23(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 23);
                            prev_cnt_v_addr25 = cnt;
                            prev_t_Start_v_addr25 = t_Start;
                            prev_t_End_v_addr25 = t;
                            prev_i_Start_v_addr25 = i_Start;
                            prev_i_End_v_addr25 = i;
                            prev_j_Start_v_addr25 = j_Start;
                            prev_j_End_v_addr25 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
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
void ref_q_addr9() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr9 = -1;
    uint64_t prev_t_Start_q_addr9 = -1;
    uint64_t prev_t_End_q_addr9 = -1;
    uint64_t prev_i_Start_q_addr9 = -1;
    uint64_t prev_i_End_q_addr9 = -1;
    uint64_t prev_j_Start_q_addr9 = -1;
    uint64_t prev_j_End_q_addr9 = -1;
    uint64_t prev_cnt_q_addr11 = -1;
    uint64_t prev_t_Start_q_addr11 = -1;
    uint64_t prev_t_End_q_addr11 = -1;
    uint64_t prev_i_Start_q_addr11 = -1;
    uint64_t prev_i_End_q_addr11 = -1;
    uint64_t prev_j_Start_q_addr11 = -1;
    uint64_t prev_j_End_q_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr9 != -1) {
            if ( calAddrq_addr9( t_Start - prev_t_Start_q_addr9 + prev_t_End_q_addr9, i_Start - prev_i_Start_q_addr9 + prev_i_End_q_addr9, j_Start - prev_j_Start_q_addr9 + prev_j_End_q_addr9) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr9, 9);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr11 != -1) {
            if ( calAddrq_addr11( t_Start - prev_t_Start_q_addr11 + prev_t_End_q_addr11, i_Start - prev_i_Start_q_addr11 + prev_i_End_q_addr11, j_Start - prev_j_Start_q_addr11 + prev_j_End_q_addr11) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr11, 11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 9);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 9);
                            prev_cnt_q_addr9 = cnt;
                            prev_t_Start_q_addr9 = t_Start;
                            prev_t_End_q_addr9 = t;
                            prev_i_Start_q_addr9 = i_Start;
                            prev_i_End_q_addr9 = i;
                            prev_j_Start_q_addr9 = j_Start;
                            prev_j_End_q_addr9 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 9);
                            prev_cnt_q_addr11 = cnt;
                            prev_t_Start_q_addr11 = t_Start;
                            prev_t_End_q_addr11 = t;
                            prev_i_Start_q_addr11 = i_Start;
                            prev_i_End_q_addr11 = i;
                            prev_j_Start_q_addr11 = j_Start;
                            prev_j_End_q_addr11 = j;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 9);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 9);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 9);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 9);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 9);
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
void ref_p_addr10() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr4 = -1;
    uint64_t prev_t_Start_p_addr4 = -1;
    uint64_t prev_t_End_p_addr4 = -1;
    uint64_t prev_i_Start_p_addr4 = -1;
    uint64_t prev_i_End_p_addr4 = -1;
    uint64_t prev_j_Start_p_addr4 = -1;
    uint64_t prev_j_End_p_addr4 = -1;
    uint64_t prev_cnt_p_addr5 = -1;
    uint64_t prev_t_Start_p_addr5 = -1;
    uint64_t prev_t_End_p_addr5 = -1;
    uint64_t prev_i_Start_p_addr5 = -1;
    uint64_t prev_i_End_p_addr5 = -1;
    uint64_t prev_j_Start_p_addr5 = -1;
    uint64_t prev_j_End_p_addr5 = -1;
    uint64_t prev_cnt_p_addr10 = -1;
    uint64_t prev_t_Start_p_addr10 = -1;
    uint64_t prev_t_End_p_addr10 = -1;
    uint64_t prev_i_Start_p_addr10 = -1;
    uint64_t prev_i_End_p_addr10 = -1;
    uint64_t prev_j_Start_p_addr10 = -1;
    uint64_t prev_j_End_p_addr10 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr4 != -1) {
            if ( calAddrp_addr4( t_Start - prev_t_Start_p_addr4 + prev_t_End_p_addr4, i_Start - prev_i_Start_p_addr4 + prev_i_End_p_addr4, j_Start - prev_j_Start_p_addr4 + prev_j_End_p_addr4) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr4, 4);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr5 != -1) {
            if ( calAddrp_addr5( t_Start - prev_t_Start_p_addr5 + prev_t_End_p_addr5, i_Start - prev_i_Start_p_addr5 + prev_i_End_p_addr5, j_Start - prev_j_Start_p_addr5 + prev_j_End_p_addr5) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr5, 5);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr10 != -1) {
            if ( calAddrp_addr10( t_Start - prev_t_Start_p_addr10 + prev_t_End_p_addr10, i_Start - prev_i_Start_p_addr10 + prev_i_End_p_addr10, j_Start - prev_j_Start_p_addr10 + prev_j_End_p_addr10) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr10, 10);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 10);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 10);
                            prev_cnt_p_addr4 = cnt;
                            prev_t_Start_p_addr4 = t_Start;
                            prev_t_End_p_addr4 = t;
                            prev_i_Start_p_addr4 = i_Start;
                            prev_i_End_p_addr4 = i;
                            prev_j_Start_p_addr4 = j_Start;
                            prev_j_End_p_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 10);
                            prev_cnt_p_addr5 = cnt;
                            prev_t_Start_p_addr5 = t_Start;
                            prev_t_End_p_addr5 = t;
                            prev_i_Start_p_addr5 = i_Start;
                            prev_i_End_p_addr5 = i;
                            prev_j_Start_p_addr5 = j_Start;
                            prev_j_End_p_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 10);
                            prev_cnt_p_addr10 = cnt;
                            prev_t_Start_p_addr10 = t_Start;
                            prev_t_End_p_addr10 = t;
                            prev_i_Start_p_addr10 = i_Start;
                            prev_i_End_p_addr10 = i;
                            prev_j_Start_p_addr10 = j_Start;
                            prev_j_End_p_addr10 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 10);
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
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 10);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 10);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 10);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 10);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr10(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 10);
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
void ref_q_addr11() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr9 = -1;
    uint64_t prev_t_Start_q_addr9 = -1;
    uint64_t prev_t_End_q_addr9 = -1;
    uint64_t prev_i_Start_q_addr9 = -1;
    uint64_t prev_i_End_q_addr9 = -1;
    uint64_t prev_j_Start_q_addr9 = -1;
    uint64_t prev_j_End_q_addr9 = -1;
    uint64_t prev_cnt_q_addr11 = -1;
    uint64_t prev_t_Start_q_addr11 = -1;
    uint64_t prev_t_End_q_addr11 = -1;
    uint64_t prev_i_Start_q_addr11 = -1;
    uint64_t prev_i_End_q_addr11 = -1;
    uint64_t prev_j_Start_q_addr11 = -1;
    uint64_t prev_j_End_q_addr11 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr9 != -1) {
            if ( calAddrq_addr9( t_Start - prev_t_Start_q_addr9 + prev_t_End_q_addr9, i_Start - prev_i_Start_q_addr9 + prev_i_End_q_addr9, j_Start - prev_j_Start_q_addr9 + prev_j_End_q_addr9) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr9, 9);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr11 != -1) {
            if ( calAddrq_addr11( t_Start - prev_t_Start_q_addr11 + prev_t_End_q_addr11, i_Start - prev_i_Start_q_addr11 + prev_i_End_q_addr11, j_Start - prev_j_Start_q_addr11 + prev_j_End_q_addr11) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr11, 11);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 11);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 11);
                            prev_cnt_q_addr9 = cnt;
                            prev_t_Start_q_addr9 = t_Start;
                            prev_t_End_q_addr9 = t;
                            prev_i_Start_q_addr9 = i_Start;
                            prev_i_End_q_addr9 = i;
                            prev_j_Start_q_addr9 = j_Start;
                            prev_j_End_q_addr9 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 11);
                            prev_cnt_q_addr11 = cnt;
                            prev_t_Start_q_addr11 = t_Start;
                            prev_t_End_q_addr11 = t;
                            prev_i_Start_q_addr11 = i_Start;
                            prev_i_End_q_addr11 = i;
                            prev_j_Start_q_addr11 = j_Start;
                            prev_j_End_q_addr11 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 11);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 11);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 11);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 11);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr11(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 11);
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
void ref_u_addr31() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr31 = -1;
    uint64_t prev_t_Start_u_addr31 = -1;
    uint64_t prev_t_End_u_addr31 = -1;
    uint64_t prev_i_Start_u_addr31 = -1;
    uint64_t prev_i_End_u_addr31 = -1;
    uint64_t prev_j_Start_u_addr31 = -1;
    uint64_t prev_j_End_u_addr31 = -1;
    uint64_t prev_cnt_u_addr33 = -1;
    uint64_t prev_t_Start_u_addr33 = -1;
    uint64_t prev_t_End_u_addr33 = -1;
    uint64_t prev_i_Start_u_addr33 = -1;
    uint64_t prev_i_End_u_addr33 = -1;
    uint64_t prev_j_Start_u_addr33 = -1;
    uint64_t prev_j_End_u_addr33 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr31 != -1) {
            if ( calAddru_addr31( t_Start - prev_t_Start_u_addr31 + prev_t_End_u_addr31, i_Start - prev_i_Start_u_addr31 + prev_i_End_u_addr31, j_Start - prev_j_Start_u_addr31 + prev_j_End_u_addr31) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr31, 31);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr33 != -1) {
            if ( calAddru_addr33( t_Start - prev_t_Start_u_addr33 + prev_t_End_u_addr33, i_Start - prev_i_Start_u_addr33 + prev_i_End_u_addr33, j_Start - prev_j_Start_u_addr33 + prev_j_End_u_addr33) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr33, 33);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 31);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 31);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 31);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 31);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 31);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr29( t, i) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 31);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 31);
                            prev_cnt_u_addr31 = cnt;
                            prev_t_Start_u_addr31 = t_Start;
                            prev_t_End_u_addr31 = t;
                            prev_i_Start_u_addr31 = i_Start;
                            prev_i_End_u_addr31 = i;
                            prev_j_Start_u_addr31 = j_Start;
                            prev_j_End_u_addr31 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr31(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 31);
                            prev_cnt_u_addr33 = cnt;
                            prev_t_Start_u_addr33 = t_Start;
                            prev_t_End_u_addr33 = t;
                            prev_i_Start_u_addr33 = i_Start;
                            prev_i_End_u_addr33 = i;
                            prev_j_Start_u_addr33 = j_Start;
                            prev_j_End_u_addr33 = j;
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
void ref_q_addr32() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr32 = -1;
    uint64_t prev_t_Start_q_addr32 = -1;
    uint64_t prev_t_End_q_addr32 = -1;
    uint64_t prev_i_Start_q_addr32 = -1;
    uint64_t prev_i_End_q_addr32 = -1;
    uint64_t prev_j_Start_q_addr32 = -1;
    uint64_t prev_j_End_q_addr32 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr32 != -1) {
            if ( calAddrq_addr32( t_Start - prev_t_Start_q_addr32 + prev_t_End_q_addr32, i_Start - prev_i_Start_q_addr32 + prev_i_End_q_addr32, j_Start - prev_j_Start_q_addr32 + prev_j_End_q_addr32) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr32, 32);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 32);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 32);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 32);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 32);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 32);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 32);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 32);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr32(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 32);
                            prev_cnt_q_addr32 = cnt;
                            prev_t_Start_q_addr32 = t_Start;
                            prev_t_End_q_addr32 = t;
                            prev_i_Start_q_addr32 = i_Start;
                            prev_i_End_q_addr32 = i;
                            prev_j_Start_q_addr32 = j_Start;
                            prev_j_End_q_addr32 = j;
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
void ref_u_addr33() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr31 = -1;
    uint64_t prev_t_Start_u_addr31 = -1;
    uint64_t prev_t_End_u_addr31 = -1;
    uint64_t prev_i_Start_u_addr31 = -1;
    uint64_t prev_i_End_u_addr31 = -1;
    uint64_t prev_j_Start_u_addr31 = -1;
    uint64_t prev_j_End_u_addr31 = -1;
    uint64_t prev_cnt_u_addr33 = -1;
    uint64_t prev_t_Start_u_addr33 = -1;
    uint64_t prev_t_End_u_addr33 = -1;
    uint64_t prev_i_Start_u_addr33 = -1;
    uint64_t prev_i_End_u_addr33 = -1;
    uint64_t prev_j_Start_u_addr33 = -1;
    uint64_t prev_j_End_u_addr33 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr31 != -1) {
            if ( calAddru_addr31( t_Start - prev_t_Start_u_addr31 + prev_t_End_u_addr31, i_Start - prev_i_Start_u_addr31 + prev_i_End_u_addr31, j_Start - prev_j_Start_u_addr31 + prev_j_End_u_addr31) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr31, 31);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr33 != -1) {
            if ( calAddru_addr33( t_Start - prev_t_Start_u_addr33 + prev_t_End_u_addr33, i_Start - prev_i_Start_u_addr33 + prev_i_End_u_addr33, j_Start - prev_j_Start_u_addr33 + prev_j_End_u_addr33) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr33, 33);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 33);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 33);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 33);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 33);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 33);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr29( t, i) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 33);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 33);
                            prev_cnt_u_addr31 = cnt;
                            prev_t_Start_u_addr31 = t_Start;
                            prev_t_End_u_addr31 = t;
                            prev_i_Start_u_addr31 = i_Start;
                            prev_i_End_u_addr31 = i;
                            prev_j_Start_u_addr31 = j_Start;
                            prev_j_End_u_addr31 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr33(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 33);
                            prev_cnt_u_addr33 = cnt;
                            prev_t_Start_u_addr33 = t_Start;
                            prev_t_End_u_addr33 = t;
                            prev_i_Start_u_addr33 = i_Start;
                            prev_i_End_u_addr33 = i;
                            prev_j_Start_u_addr33 = j_Start;
                            prev_j_End_u_addr33 = j;
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
void ref_v_addr16() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr14 = -1;
    uint64_t prev_t_Start_v_addr14 = -1;
    uint64_t prev_t_End_v_addr14 = -1;
    uint64_t prev_i_Start_v_addr14 = -1;
    uint64_t prev_i_End_v_addr14 = -1;
    uint64_t prev_j_Start_v_addr14 = -1;
    uint64_t prev_j_End_v_addr14 = -1;
    uint64_t prev_cnt_v_addr16 = -1;
    uint64_t prev_t_Start_v_addr16 = -1;
    uint64_t prev_t_End_v_addr16 = -1;
    uint64_t prev_i_Start_v_addr16 = -1;
    uint64_t prev_i_End_v_addr16 = -1;
    uint64_t prev_j_Start_v_addr16 = -1;
    uint64_t prev_j_End_v_addr16 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr14 != -1) {
            if ( calAddrv_addr14( t_Start - prev_t_Start_v_addr14 + prev_t_End_v_addr14, i_Start - prev_i_Start_v_addr14 + prev_i_End_v_addr14, j_Start - prev_j_Start_v_addr14 + prev_j_End_v_addr14) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr14, 14);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr16 != -1) {
            if ( calAddrv_addr16( t_Start - prev_t_Start_v_addr16 + prev_t_End_v_addr16, i_Start - prev_i_Start_v_addr16 + prev_i_End_v_addr16, j_Start - prev_j_Start_v_addr16 + prev_j_End_v_addr16) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr16, 16);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 16);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 16);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr12( t, i) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 16);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 16);
                            prev_cnt_v_addr14 = cnt;
                            prev_t_Start_v_addr14 = t_Start;
                            prev_t_End_v_addr14 = t;
                            prev_i_Start_v_addr14 = i_Start;
                            prev_i_End_v_addr14 = i;
                            prev_j_Start_v_addr14 = j_Start;
                            prev_j_End_v_addr14 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 16);
                            prev_cnt_v_addr16 = cnt;
                            prev_t_Start_v_addr16 = t_Start;
                            prev_t_End_v_addr16 = t;
                            prev_i_Start_v_addr16 = i_Start;
                            prev_i_End_v_addr16 = i;
                            prev_j_Start_v_addr16 = j_Start;
                            prev_j_End_v_addr16 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 16);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 16);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr16(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 16);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
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
void ref_p_addr18() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr18 = -1;
    uint64_t prev_t_Start_p_addr18 = -1;
    uint64_t prev_t_End_p_addr18 = -1;
    uint64_t prev_i_Start_p_addr18 = -1;
    uint64_t prev_i_End_p_addr18 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr18 != -1) {
            if ( calAddrp_addr18( t_Start - prev_t_Start_p_addr18 + prev_t_End_p_addr18, i_Start - prev_i_Start_p_addr18 + prev_i_End_p_addr18) == calAddrp_addr18(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_p_addr18, 18);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr18(t_Start, i_Start)) {
                        rtHistoCal(cnt, 18);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                            rtHistoCal(cnt, 18);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                            rtHistoCal(cnt, 18);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                            rtHistoCal(cnt, 18);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                            rtHistoCal(cnt, 18);
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr18(t_Start, i_Start)) {
                        rtHistoCal(cnt, 18);
                        prev_cnt_p_addr18 = cnt;
                        prev_t_Start_p_addr18 = t_Start;
                        prev_t_End_p_addr18 = t;
                        prev_i_Start_p_addr18 = i_Start;
                        prev_i_End_p_addr18 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                            rtHistoCal(cnt, 18);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                            rtHistoCal(cnt, 18);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                            rtHistoCal(cnt, 18);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr18(t_Start, i_Start)) {
                            rtHistoCal(cnt, 18);
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
void ref_u_addr19() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr17 = -1;
    uint64_t prev_t_Start_u_addr17 = -1;
    uint64_t prev_t_End_u_addr17 = -1;
    uint64_t prev_i_Start_u_addr17 = -1;
    uint64_t prev_i_End_u_addr17 = -1;
    uint64_t prev_cnt_u_addr19 = -1;
    uint64_t prev_t_Start_u_addr19 = -1;
    uint64_t prev_t_End_u_addr19 = -1;
    uint64_t prev_i_Start_u_addr19 = -1;
    uint64_t prev_i_End_u_addr19 = -1;
    uint64_t prev_cnt_u_addr29 = -1;
    uint64_t prev_t_Start_u_addr29 = -1;
    uint64_t prev_t_End_u_addr29 = -1;
    uint64_t prev_i_Start_u_addr29 = -1;
    uint64_t prev_i_End_u_addr29 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr17 != -1) {
            if ( calAddru_addr17( t_Start - prev_t_Start_u_addr17 + prev_t_End_u_addr17, i_Start - prev_i_Start_u_addr17 + prev_i_End_u_addr17) == calAddru_addr19(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr17, 17);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr19 != -1) {
            if ( calAddru_addr19( t_Start - prev_t_Start_u_addr19 + prev_t_End_u_addr19, i_Start - prev_i_Start_u_addr19 + prev_i_End_u_addr19) == calAddru_addr19(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr19, 19);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr29 != -1) {
            if ( calAddru_addr29( t_Start - prev_t_Start_u_addr29 + prev_t_End_u_addr29, i_Start - prev_i_Start_u_addr29 + prev_i_End_u_addr29) == calAddru_addr19(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr29, 29);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                            rtHistoCal(cnt, 19);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                            rtHistoCal(cnt, 19);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                            rtHistoCal(cnt, 19);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr19(t_Start, i_Start)) {
                        rtHistoCal(cnt, 19);
                        prev_cnt_u_addr17 = cnt;
                        prev_t_Start_u_addr17 = t_Start;
                        prev_t_End_u_addr17 = t;
                        prev_i_Start_u_addr17 = i_Start;
                        prev_i_End_u_addr17 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr19(t_Start, i_Start)) {
                        rtHistoCal(cnt, 19);
                        prev_cnt_u_addr19 = cnt;
                        prev_t_Start_u_addr19 = t_Start;
                        prev_t_End_u_addr19 = t;
                        prev_i_Start_u_addr19 = i_Start;
                        prev_i_End_u_addr19 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr29( t, i) == calAddru_addr19(t_Start, i_Start)) {
                        rtHistoCal(cnt, 19);
                        prev_cnt_u_addr29 = cnt;
                        prev_t_Start_u_addr29 = t_Start;
                        prev_t_End_u_addr29 = t;
                        prev_i_Start_u_addr29 = i_Start;
                        prev_i_End_u_addr29 = i;
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                            rtHistoCal(cnt, 19);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr19(t_Start, i_Start)) {
                            rtHistoCal(cnt, 19);
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
void ref_v_addr24() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr23 = -1;
    uint64_t prev_t_Start_v_addr23 = -1;
    uint64_t prev_t_End_v_addr23 = -1;
    uint64_t prev_i_Start_v_addr23 = -1;
    uint64_t prev_i_End_v_addr23 = -1;
    uint64_t prev_j_Start_v_addr23 = -1;
    uint64_t prev_j_End_v_addr23 = -1;
    uint64_t prev_cnt_v_addr24 = -1;
    uint64_t prev_t_Start_v_addr24 = -1;
    uint64_t prev_t_End_v_addr24 = -1;
    uint64_t prev_i_Start_v_addr24 = -1;
    uint64_t prev_i_End_v_addr24 = -1;
    uint64_t prev_j_Start_v_addr24 = -1;
    uint64_t prev_j_End_v_addr24 = -1;
    uint64_t prev_cnt_v_addr25 = -1;
    uint64_t prev_t_Start_v_addr25 = -1;
    uint64_t prev_t_End_v_addr25 = -1;
    uint64_t prev_i_Start_v_addr25 = -1;
    uint64_t prev_i_End_v_addr25 = -1;
    uint64_t prev_j_Start_v_addr25 = -1;
    uint64_t prev_j_End_v_addr25 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr23 != -1) {
            if ( calAddrv_addr23( t_Start - prev_t_Start_v_addr23 + prev_t_End_v_addr23, i_Start - prev_i_Start_v_addr23 + prev_i_End_v_addr23, j_Start - prev_j_Start_v_addr23 + prev_j_End_v_addr23) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr23, 23);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr24 != -1) {
            if ( calAddrv_addr24( t_Start - prev_t_Start_v_addr24 + prev_t_End_v_addr24, i_Start - prev_i_Start_v_addr24 + prev_i_End_v_addr24, j_Start - prev_j_Start_v_addr24 + prev_j_End_v_addr24) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr24, 24);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr25 != -1) {
            if ( calAddrv_addr25( t_Start - prev_t_Start_v_addr25 + prev_t_End_v_addr25, i_Start - prev_i_Start_v_addr25 + prev_i_End_v_addr25, j_Start - prev_j_Start_v_addr25 + prev_j_End_v_addr25) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr25, 25);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 24);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 24);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr12( t, i) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 24);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 24);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 24);
                            goto EndSample;
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
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 24);
                            prev_cnt_v_addr23 = cnt;
                            prev_t_Start_v_addr23 = t_Start;
                            prev_t_End_v_addr23 = t;
                            prev_i_Start_v_addr23 = i_Start;
                            prev_i_End_v_addr23 = i;
                            prev_j_Start_v_addr23 = j_Start;
                            prev_j_End_v_addr23 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 24);
                            prev_cnt_v_addr24 = cnt;
                            prev_t_Start_v_addr24 = t_Start;
                            prev_t_End_v_addr24 = t;
                            prev_i_Start_v_addr24 = i_Start;
                            prev_i_End_v_addr24 = i;
                            prev_j_Start_v_addr24 = j_Start;
                            prev_j_End_v_addr24 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr24(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 24);
                            prev_cnt_v_addr25 = cnt;
                            prev_t_Start_v_addr25 = t_Start;
                            prev_t_End_v_addr25 = t;
                            prev_i_Start_v_addr25 = i_Start;
                            prev_i_End_v_addr25 = i;
                            prev_j_Start_v_addr25 = j_Start;
                            prev_j_End_v_addr25 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
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
void ref_v_addr25() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr23 = -1;
    uint64_t prev_t_Start_v_addr23 = -1;
    uint64_t prev_t_End_v_addr23 = -1;
    uint64_t prev_i_Start_v_addr23 = -1;
    uint64_t prev_i_End_v_addr23 = -1;
    uint64_t prev_j_Start_v_addr23 = -1;
    uint64_t prev_j_End_v_addr23 = -1;
    uint64_t prev_cnt_v_addr24 = -1;
    uint64_t prev_t_Start_v_addr24 = -1;
    uint64_t prev_t_End_v_addr24 = -1;
    uint64_t prev_i_Start_v_addr24 = -1;
    uint64_t prev_i_End_v_addr24 = -1;
    uint64_t prev_j_Start_v_addr24 = -1;
    uint64_t prev_j_End_v_addr24 = -1;
    uint64_t prev_cnt_v_addr25 = -1;
    uint64_t prev_t_Start_v_addr25 = -1;
    uint64_t prev_t_End_v_addr25 = -1;
    uint64_t prev_i_Start_v_addr25 = -1;
    uint64_t prev_i_End_v_addr25 = -1;
    uint64_t prev_j_Start_v_addr25 = -1;
    uint64_t prev_j_End_v_addr25 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr23 != -1) {
            if ( calAddrv_addr23( t_Start - prev_t_Start_v_addr23 + prev_t_End_v_addr23, i_Start - prev_i_Start_v_addr23 + prev_i_End_v_addr23, j_Start - prev_j_Start_v_addr23 + prev_j_End_v_addr23) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr23, 23);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr24 != -1) {
            if ( calAddrv_addr24( t_Start - prev_t_Start_v_addr24 + prev_t_End_v_addr24, i_Start - prev_i_Start_v_addr24 + prev_i_End_v_addr24, j_Start - prev_j_Start_v_addr24 + prev_j_End_v_addr24) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr24, 24);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr25 != -1) {
            if ( calAddrv_addr25( t_Start - prev_t_Start_v_addr25 + prev_t_End_v_addr25, i_Start - prev_i_Start_v_addr25 + prev_i_End_v_addr25, j_Start - prev_j_Start_v_addr25 + prev_j_End_v_addr25) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr25, 25);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 25);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 25);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr12( t, i) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 25);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr14( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 25);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr16( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 25);
                            goto EndSample;
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
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr23( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 25);
                            prev_cnt_v_addr23 = cnt;
                            prev_t_Start_v_addr23 = t_Start;
                            prev_t_End_v_addr23 = t;
                            prev_i_Start_v_addr23 = i_Start;
                            prev_i_End_v_addr23 = i;
                            prev_j_Start_v_addr23 = j_Start;
                            prev_j_End_v_addr23 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr24( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 25);
                            prev_cnt_v_addr24 = cnt;
                            prev_t_Start_v_addr24 = t_Start;
                            prev_t_End_v_addr24 = t;
                            prev_i_Start_v_addr24 = i_Start;
                            prev_i_End_v_addr24 = i;
                            prev_j_Start_v_addr24 = j_Start;
                            prev_j_End_v_addr24 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr25( t, i, j) == calAddrv_addr25(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 25);
                            prev_cnt_v_addr25 = cnt;
                            prev_t_Start_v_addr25 = t_Start;
                            prev_t_End_v_addr25 = t;
                            prev_i_Start_v_addr25 = i_Start;
                            prev_i_End_v_addr25 = i;
                            prev_j_Start_v_addr25 = j_Start;
                            prev_j_End_v_addr25 = j;
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
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
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
void ref_q_addr26() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr26 = -1;
    uint64_t prev_t_Start_q_addr26 = -1;
    uint64_t prev_t_End_q_addr26 = -1;
    uint64_t prev_i_Start_q_addr26 = -1;
    uint64_t prev_i_End_q_addr26 = -1;
    uint64_t prev_j_Start_q_addr26 = -1;
    uint64_t prev_j_End_q_addr26 = -1;
    uint64_t prev_cnt_q_addr28 = -1;
    uint64_t prev_t_Start_q_addr28 = -1;
    uint64_t prev_t_End_q_addr28 = -1;
    uint64_t prev_i_Start_q_addr28 = -1;
    uint64_t prev_i_End_q_addr28 = -1;
    uint64_t prev_j_Start_q_addr28 = -1;
    uint64_t prev_j_End_q_addr28 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr26 != -1) {
            if ( calAddrq_addr26( t_Start - prev_t_Start_q_addr26 + prev_t_End_q_addr26, i_Start - prev_i_Start_q_addr26 + prev_i_End_q_addr26, j_Start - prev_j_Start_q_addr26 + prev_j_End_q_addr26) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr26, 26);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr28 != -1) {
            if ( calAddrq_addr28( t_Start - prev_t_Start_q_addr28 + prev_t_End_q_addr28, i_Start - prev_i_Start_q_addr28 + prev_i_End_q_addr28, j_Start - prev_j_Start_q_addr28 + prev_j_End_q_addr28) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr28, 28);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 26);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 26);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 26);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 26);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 26);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 26);
                            prev_cnt_q_addr26 = cnt;
                            prev_t_Start_q_addr26 = t_Start;
                            prev_t_End_q_addr26 = t;
                            prev_i_Start_q_addr26 = i_Start;
                            prev_i_End_q_addr26 = i;
                            prev_j_Start_q_addr26 = j_Start;
                            prev_j_End_q_addr26 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 26);
                            prev_cnt_q_addr28 = cnt;
                            prev_t_Start_q_addr28 = t_Start;
                            prev_t_End_q_addr28 = t;
                            prev_i_Start_q_addr28 = i_Start;
                            prev_i_End_q_addr28 = i;
                            prev_j_Start_q_addr28 = j_Start;
                            prev_j_End_q_addr28 = j;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr26(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 26);
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
void ref_p_addr27() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr21 = -1;
    uint64_t prev_t_Start_p_addr21 = -1;
    uint64_t prev_t_End_p_addr21 = -1;
    uint64_t prev_i_Start_p_addr21 = -1;
    uint64_t prev_i_End_p_addr21 = -1;
    uint64_t prev_j_Start_p_addr21 = -1;
    uint64_t prev_j_End_p_addr21 = -1;
    uint64_t prev_cnt_p_addr22 = -1;
    uint64_t prev_t_Start_p_addr22 = -1;
    uint64_t prev_t_End_p_addr22 = -1;
    uint64_t prev_i_Start_p_addr22 = -1;
    uint64_t prev_i_End_p_addr22 = -1;
    uint64_t prev_j_Start_p_addr22 = -1;
    uint64_t prev_j_End_p_addr22 = -1;
    uint64_t prev_cnt_p_addr27 = -1;
    uint64_t prev_t_Start_p_addr27 = -1;
    uint64_t prev_t_End_p_addr27 = -1;
    uint64_t prev_i_Start_p_addr27 = -1;
    uint64_t prev_i_End_p_addr27 = -1;
    uint64_t prev_j_Start_p_addr27 = -1;
    uint64_t prev_j_End_p_addr27 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr21 != -1) {
            if ( calAddrp_addr21( t_Start - prev_t_Start_p_addr21 + prev_t_End_p_addr21, i_Start - prev_i_Start_p_addr21 + prev_i_End_p_addr21, j_Start - prev_j_Start_p_addr21 + prev_j_End_p_addr21) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr21, 21);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr22 != -1) {
            if ( calAddrp_addr22( t_Start - prev_t_Start_p_addr22 + prev_t_End_p_addr22, i_Start - prev_i_Start_p_addr22 + prev_i_End_p_addr22, j_Start - prev_j_Start_p_addr22 + prev_j_End_p_addr22) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr22, 22);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr27 != -1) {
            if ( calAddrp_addr27( t_Start - prev_t_Start_p_addr27 + prev_t_End_p_addr27, i_Start - prev_i_Start_p_addr27 + prev_i_End_p_addr27, j_Start - prev_j_Start_p_addr27 + prev_j_End_p_addr27) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr27, 27);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 27);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 27);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 27);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 27);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 27);
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 27);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 27);
                            prev_cnt_p_addr21 = cnt;
                            prev_t_Start_p_addr21 = t_Start;
                            prev_t_End_p_addr21 = t;
                            prev_i_Start_p_addr21 = i_Start;
                            prev_i_End_p_addr21 = i;
                            prev_j_Start_p_addr21 = j_Start;
                            prev_j_End_p_addr21 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 27);
                            prev_cnt_p_addr22 = cnt;
                            prev_t_Start_p_addr22 = t_Start;
                            prev_t_End_p_addr22 = t;
                            prev_i_Start_p_addr22 = i_Start;
                            prev_i_End_p_addr22 = i;
                            prev_j_Start_p_addr22 = j_Start;
                            prev_j_End_p_addr22 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 27);
                            prev_cnt_p_addr27 = cnt;
                            prev_t_Start_p_addr27 = t_Start;
                            prev_t_End_p_addr27 = t;
                            prev_i_Start_p_addr27 = i_Start;
                            prev_i_End_p_addr27 = i;
                            prev_j_Start_p_addr27 = j_Start;
                            prev_j_End_p_addr27 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr27(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 27);
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
void ref_q_addr28() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr26 = -1;
    uint64_t prev_t_Start_q_addr26 = -1;
    uint64_t prev_t_End_q_addr26 = -1;
    uint64_t prev_i_Start_q_addr26 = -1;
    uint64_t prev_i_End_q_addr26 = -1;
    uint64_t prev_j_Start_q_addr26 = -1;
    uint64_t prev_j_End_q_addr26 = -1;
    uint64_t prev_cnt_q_addr28 = -1;
    uint64_t prev_t_Start_q_addr28 = -1;
    uint64_t prev_t_End_q_addr28 = -1;
    uint64_t prev_i_Start_q_addr28 = -1;
    uint64_t prev_i_End_q_addr28 = -1;
    uint64_t prev_j_Start_q_addr28 = -1;
    uint64_t prev_j_End_q_addr28 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr26 != -1) {
            if ( calAddrq_addr26( t_Start - prev_t_Start_q_addr26 + prev_t_End_q_addr26, i_Start - prev_i_Start_q_addr26 + prev_i_End_q_addr26, j_Start - prev_j_Start_q_addr26 + prev_j_End_q_addr26) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr26, 26);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr28 != -1) {
            if ( calAddrq_addr28( t_Start - prev_t_Start_q_addr28 + prev_t_End_q_addr28, i_Start - prev_i_Start_q_addr28 + prev_i_End_q_addr28, j_Start - prev_j_Start_q_addr28 + prev_j_End_q_addr28) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr28, 28);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr3( t, i) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 28);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr9( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 28);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr11( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 28);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr15( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 28);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr20( t, i) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 28);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr26( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 28);
                            prev_cnt_q_addr26 = cnt;
                            prev_t_Start_q_addr26 = t_Start;
                            prev_t_End_q_addr26 = t;
                            prev_i_Start_q_addr26 = i_Start;
                            prev_i_End_q_addr26 = i;
                            prev_j_Start_q_addr26 = j_Start;
                            prev_j_End_q_addr26 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr28( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 28);
                            prev_cnt_q_addr28 = cnt;
                            prev_t_Start_q_addr28 = t_Start;
                            prev_t_End_q_addr28 = t;
                            prev_i_Start_q_addr28 = i_Start;
                            prev_i_End_q_addr28 = i;
                            prev_j_Start_q_addr28 = j_Start;
                            prev_j_End_q_addr28 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr32( t, i, j) == calAddrq_addr28(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 28);
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
void ref_u_addr29() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr17 = -1;
    uint64_t prev_t_Start_u_addr17 = -1;
    uint64_t prev_t_End_u_addr17 = -1;
    uint64_t prev_i_Start_u_addr17 = -1;
    uint64_t prev_i_End_u_addr17 = -1;
    uint64_t prev_cnt_u_addr19 = -1;
    uint64_t prev_t_Start_u_addr19 = -1;
    uint64_t prev_t_End_u_addr19 = -1;
    uint64_t prev_i_Start_u_addr19 = -1;
    uint64_t prev_i_End_u_addr19 = -1;
    uint64_t prev_cnt_u_addr29 = -1;
    uint64_t prev_t_Start_u_addr29 = -1;
    uint64_t prev_t_End_u_addr29 = -1;
    uint64_t prev_i_Start_u_addr29 = -1;
    uint64_t prev_i_End_u_addr29 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 25;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr17 != -1) {
            if ( calAddru_addr17( t_Start - prev_t_Start_u_addr17 + prev_t_End_u_addr17, i_Start - prev_i_Start_u_addr17 + prev_i_End_u_addr17) == calAddru_addr29(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr17, 17);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr19 != -1) {
            if ( calAddru_addr19( t_Start - prev_t_Start_u_addr19 + prev_t_End_u_addr19, i_Start - prev_i_Start_u_addr19 + prev_i_End_u_addr19) == calAddru_addr29(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr19, 19);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr29 != -1) {
            if ( calAddru_addr29( t_Start - prev_t_Start_u_addr29 + prev_t_End_u_addr29, i_Start - prev_i_Start_u_addr29 + prev_i_End_u_addr29) == calAddru_addr29(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr29, 29);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                            rtHistoCal(cnt, 29);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                            rtHistoCal(cnt, 29);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr8( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                            rtHistoCal(cnt, 29);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr17( t, i) == calAddru_addr29(t_Start, i_Start)) {
                        rtHistoCal(cnt, 29);
                        prev_cnt_u_addr17 = cnt;
                        prev_t_Start_u_addr17 = t_Start;
                        prev_t_End_u_addr17 = t;
                        prev_i_Start_u_addr17 = i_Start;
                        prev_i_End_u_addr17 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr19( t, i) == calAddru_addr29(t_Start, i_Start)) {
                        rtHistoCal(cnt, 29);
                        prev_cnt_u_addr19 = cnt;
                        prev_t_Start_u_addr19 = t_Start;
                        prev_t_End_u_addr19 = t;
                        prev_i_Start_u_addr19 = i_Start;
                        prev_i_End_u_addr19 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr29( t, i) == calAddru_addr29(t_Start, i_Start)) {
                        rtHistoCal(cnt, 29);
                        prev_cnt_u_addr29 = cnt;
                        prev_t_Start_u_addr29 = t_Start;
                        prev_t_End_u_addr29 = t;
                        prev_i_Start_u_addr29 = i_Start;
                        prev_i_End_u_addr29 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr31( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                            rtHistoCal(cnt, 29);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr33( t, i, j) == calAddru_addr29(t_Start, i_Start)) {
                            rtHistoCal(cnt, 29);
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
void ref_p_addr30() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr30 = -1;
    uint64_t prev_t_Start_p_addr30 = -1;
    uint64_t prev_t_End_p_addr30 = -1;
    uint64_t prev_i_Start_p_addr30 = -1;
    uint64_t prev_i_End_p_addr30 = -1;
    uint64_t prev_j_Start_p_addr30 = -1;
    uint64_t prev_j_End_p_addr30 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1305;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1022 - 1 + 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr30 != -1) {
            if ( calAddrp_addr30( t_Start - prev_t_Start_p_addr30 + prev_t_End_p_addr30, i_Start - prev_i_Start_p_addr30 + prev_i_End_p_addr30, j_Start - prev_j_Start_p_addr30 + prev_j_End_p_addr30) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr30, 30);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr1( t, i) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 30);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 30);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr5( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 30);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr10( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 30);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr13( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 30);
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
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr18( t, i) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt, 30);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr21( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 30);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr22( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 30);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr27( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 30);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr30( t, i, j) == calAddrp_addr30(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt, 30);
                            prev_cnt_p_addr30 = cnt;
                            prev_t_Start_p_addr30 = t_Start;
                            prev_t_End_p_addr30 = t;
                            prev_i_Start_p_addr30 = i_Start;
                            prev_i_End_p_addr30 = i;
                            prev_j_Start_p_addr30 = j_Start;
                            prev_j_End_p_addr30 = j;
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
    ref_p_addr4();
    ref_u_addr7();
    ref_u_addr8();
    ref_v_addr0();
    ref_p_addr1();
    ref_v_addr2();
    ref_q_addr3();
    ref_v_addr12();
    ref_p_addr13();
    ref_v_addr14();
    ref_q_addr15();
    ref_p_addr5();
    ref_u_addr17();
    ref_u_addr6();
    ref_q_addr20();
    ref_p_addr21();
    ref_p_addr22();
    ref_v_addr23();
    ref_q_addr9();
    ref_p_addr10();
    ref_q_addr11();
    ref_u_addr31();
    ref_q_addr32();
    ref_u_addr33();
    ref_v_addr16();
    ref_p_addr18();
    ref_u_addr19();
    ref_v_addr24();
    ref_v_addr25();
    ref_q_addr26();
    ref_p_addr27();
    ref_q_addr28();
    ref_u_addr29();
    ref_p_addr30();
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
 /* Analyze function: adi */ 
