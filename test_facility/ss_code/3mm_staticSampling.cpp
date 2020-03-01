
 /* Start to analysis array index
Array index info: Total number of references: 15
E.addr ((i * 128) + j)
E.addr ((i * 128) + j)
E.addr ((i * 128) + j)
C.addr ((i * 128) + k)
D.addr ((k * 128) + j)
A.addr ((i * 128) + k)
B.addr ((k * 128) + j)
G.addr ((i * 128) + j)
E.addr ((i * 128) + k)
F.addr ((k * 128) + j)
F.addr ((i * 128) + j)
F.addr ((i * 128) + j)
F.addr ((i * 128) + j)
G.addr ((i * 128) + j)
G.addr ((i * 128) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %ni
i32 %nj
i32 %nk
i32 %nl
i32 %nm
double* %E
double* %A
double* %B
double* %F
double* %C
double* %D
double* %G

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 128)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 128)
----Loop inc: (j + 1)
----Loop predicate: <
------array access E.addr ((i * 128) + j)
------k
------Loop Bound: (0, 128)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 128) + k)
--------array access B.addr ((k * 128) + j)
--------array access E.addr ((i * 128) + j)
--------array access E.addr ((i * 128) + j)
--i
--Loop Bound: (0, 128)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 128)
----Loop inc: (j + 1)
----Loop predicate: <
------array access F.addr ((i * 128) + j)
------k
------Loop Bound: (0, 128)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access C.addr ((i * 128) + k)
--------array access D.addr ((k * 128) + j)
--------array access F.addr ((i * 128) + j)
--------array access F.addr ((i * 128) + j)
--i
--Loop Bound: (0, 128)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 128)
----Loop inc: (j + 1)
----Loop predicate: <
------array access G.addr ((i * 128) + j)
------k
------Loop Bound: (0, 128)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access E.addr ((i * 128) + k)
--------array access F.addr ((k * 128) + j)
--------array access G.addr ((i * 128) + j)
--------array access G.addr ((i * 128) + j)

Finish analysis loops */ 
/* # of Out-most Loops: 3 */ 

 /* Start transform loop tree
	for.cond
	for.cond1
	for.cond4
	for.cond27
	for.cond30
	for.cond37
	for.cond63
	for.cond66
	for.cond73
----------------
--|  LoopNode  |
----------------
------------------
----|  LoopNode  |
------------------
--------------------
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------
--------------------
------|  LoopNode  |
--------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------
--|  LoopNode  |
----------------
------------------
----|  LoopNode  |
------------------
--------------------
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------
--------------------
------|  LoopNode  |
--------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------
--|  LoopNode  |
----------------
------------------
----|  LoopNode  |
------------------
--------------------
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------
--------------------
------|  LoopNode  |
--------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------

Finish transform loop tree */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 12
------Sample number: 163
--------Sample number: 2097
----Sample number: 12
------Sample number: 163
--------Sample number: 2097
----Sample number: 12
------Sample number: 163
--------Sample number: 2097
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
#include <map>
#include <set>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <numeric>
#include <utility>
#include <thread>
#include <mutex>
#include <algorithm>
#ifndef THREAD_NUM
#    define THREAD_NUM   4
#endif
#ifndef BIN_SIZE
#    define BIN_SIZE   4
#endif
#ifndef CHUNK_SIZE
#    define CHUNK_SIZE   4
#endif
using namespace std;
std::mutex mtx;
std::map<string, map<uint64_t, double>> refRT;
std::map<uint64_t, double> RT;
std::map<uint64_t, double> MR;
void rtHistoCal( map<uint64_t, double> &rth, uint64_t rt, double val ) {
    if ( val <= 0) {
;        return;
    }
    std::unique_lock<std::mutex> lck (mtx,std::defer_lock);
    lck.lock();
    if (rth.find(rt) == rth.end()) { 
        rth[rt] = val;
    } else {
        rth[rt] += val;
    }
    lck.unlock();
    return;
}

void refRTHistoCal(map<string, map<uint64_t, double>> &rth, uint64_t rt, double val, string ref ) {
    if ( val <= 0) {
;        return;
    }
    unique_lock< mutex> lck (mtx, defer_lock);
    lck.lock();
    if (rth.find(ref) == rth.end()) { 
        rth[ref][rt] = val;
    }
    else {
        if (rth[ref].find(rt) == rth[ref].end()) { 
            rth[ref][rt] = val;
        } else {
            rth[ref][rt] += val;
        }
    }
    lck.unlock();
    return;
}
void subBlkRT(map<uint64_t, double> &rth, int rt, double cnt) {
    int msb = 0;
    int tmp_rt = rt;
    while(tmp_rt != 0) {
        tmp_rt = tmp_rt / 2;
        ++msb;
    }
    if (msb >= BIN_SIZE) {
        int diff = (pow(2, msb) - pow(2, msb-1)) / BIN_SIZE;
        for (int b = pow(2, msb-1); b <= pow(2, msb); b+=diff) {
            if (rt < b) {
                rtHistoCal(rth, b - diff, cnt);
                break;
            }
        }
    }
    else {
        rtHistoCal(rth, pow(2, msb-1), cnt);
    }
    return;
}

void refSubBlkRT(map<string, map<uint64_t, double>> &rth, uint64_t rt, double cnt, string ref) {
    int msb = 0;
    int tmp_rt = rt;
    while(tmp_rt != 0) {
        tmp_rt = tmp_rt / 2;
        ++msb;
    }
    if (msb >= BIN_SIZE) {
        int diff = (pow(2, msb) - pow(2, msb-1)) / BIN_SIZE;
        for (int b = pow(2, msb-1); b <= pow(2, msb); b+=diff) {
            if (rt < b) {
                refRTHistoCal(rth, b - diff, cnt, ref);
                break;
            }
        }
    }
    else {
        refRTHistoCal(rth, pow(2, msb-1), cnt, ref);
    }
    return;
}
/* Merge the refRT to RT */
void rtMerge() {
    for(map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {
        for (map<uint64_t, double>::iterator iit = it->second.begin(); iit != it->second.end(); ++iit) {
            rtHistoCal(RT, iit->first, iit->second);
        }
    }
    return;
}
/* Smoothing the per-reference RTHisto based on Gasussian */
void gaussian_smoothing(map<uint64_t, double> &rth, bool scale, double sigma) {
    map<uint64_t, double> tmp;
    double sum = accumulate(begin(rth), end(rth), 0.0, [](const double previous, const pair<uint64_t, double>& p) { return previous + p.second; });
    double sum_P = 0.0;
    for(map<uint64_t, double>::iterator it = rth.begin(); it != rth.end(); ++it) {
        sum_P += it->second;
        uint64_t mu = it->first;
        if (sum_P >= 0.9 * sum) {
            rtHistoCal(tmp, mu, it->second);
            continue;
        }
        if (scale) { mu = mu * THREAD_NUM; }
        uint64_t start_b = ( mu / THREAD_NUM) >= 1 ? mu / THREAD_NUM : 1;
        uint64_t end_b = mu * THREAD_NUM;
        for(uint64_t b = start_b; b <= end_b; b++) {
            double c =  (1 /  (sqrt(2 * M_PI) * sigma));
            double val = it->second * c * exp( -1 * pow((b - mu), 2.0) / (2 * pow(sigma, 2.0)));
            rtHistoCal(tmp, b, val);
        }
    }
    rth.clear();
    for(map<uint64_t, double>::iterator it = tmp.begin(); it != tmp.end(); ++it) {
        subBlkRT(rth, it->first, it->second);
    }
    return;
}
/* Smoothing the per-reference RTHisto Uniformly. Equally splite the RT to a range of bins */
void uniform_smoothing(map<uint64_t, double> &rth, bool scale) {
    map<uint64_t, double> tmp;
    double sum = accumulate(begin(rth), end(rth), 0.0, [](const double previous, const pair<uint64_t, double>& p) { return previous + p.second; });
    double sum_P = 0.0;
    for(map<uint64_t, double>::iterator it = rth.begin(); it != rth.end(); ++it) {
        uint64_t mu = it->first;
        sum_P += it->second;
        /* Do uniform distribution for all ri, distribute from ri / THREAD_NUM to ri * THREAD_NUM  */
        if (it->second <= 0.0) { continue; }
        if (sum_P >= 0.9 * sum ) { 
            rtHistoCal(tmp, mu, it->second);
            continue;
        }
        if (scale) { mu = mu * THREAD_NUM; }
        uint64_t start = (mu / THREAD_NUM) >= 1 ? mu / THREAD_NUM : 1;
        uint64_t end = mu * THREAD_NUM;
        double split_val = it->second / (end - start + 1);
        for (int b = start; b <= end; b++) {
            rtHistoCal(tmp, b, split_val);
        }
    }
    rth.clear();
    for(map<uint64_t, double>::iterator it = tmp.begin(); it != tmp.end(); ++it) {
        subBlkRT(rth, it->first, it->second);
    }
    return;
}
void group_gaussian_smoothing(double sigma) {
    vector<string> vec = { "", "" };
    for (map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {
        gaussian_smoothing(it->second, find(vec.begin(), vec.end(), it->first) != vec.end(), sigma);
    }
}
/* Array E_addr	i j */ 
/* Array E_addr	i j */ 
/* Array E_addr	i j */ 
/* Array C_addr	i k */ 
/* Array D_addr	k j */ 
/* Array A_addr	i k */ 
/* Array B_addr	k j */ 
/* Array G_addr	i j */ 
/* Array E_addr	i k */ 
/* Array F_addr	k j */ 
/* Array F_addr	i j */ 
/* Array F_addr	i j */ 
/* Array F_addr	i j */ 
/* Array G_addr	i j */ 
/* Array G_addr	i j */ 
void group_uniform_smoothing() {
    vector<string> vec = {"D_addr0", "B_addr0", "F_addr3" };
    for (map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {
    /* Low Sensitive Program - Scaling outmost independent array references only */
        uniform_smoothing(it->second, find(vec.begin(), vec.end(), it->first) != vec.end());
    }
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
        cout << it->first << ", " << it->second << "\n";
    }
    return;
}

void refRTDump() {
    for (map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {
        cout << "Start to dump reuse time histogram for " << it->first << "\n";
        for (map<uint64_t, double>::iterator iit = it->second.begin(), eiit = it->second.end(); iit != eiit; ++iit) {
            cout << it->first << "," << iit->first << "," << iit->second << endl;
        }
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
        cout << it1->first << ", " << it1->second << endl;
        if (it1 != it2) {
            cout << it2->first << ", " << it2->second << endl;
        }
        it1 = ++it2;
        it2 = it1;
    }
    return;
}
/* E_addr ((i * 128) + j) 0 */
int calAddrE_addr0( int i, int j) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 128) + k) 0 */
int calAddrA_addr0( int i, int j, int k) {
    int result = (((i * 128) + k)) * 8 / 64;
    return result;
}
/* B_addr ((k * 128) + j) 0 */
int calAddrB_addr0( int i, int j, int k) {
    int result = (((k * 128) + j)) * 8 / 64;
    return result;
}
/* E_addr ((i * 128) + j) 1 */
int calAddrE_addr1( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* E_addr ((i * 128) + j) 2 */
int calAddrE_addr2( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* F_addr ((i * 128) + j) 0 */
int calAddrF_addr0( int i, int j) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* C_addr ((i * 128) + k) 0 */
int calAddrC_addr0( int i, int j, int k) {
    int result = (((i * 128) + k)) * 8 / 64;
    return result;
}
/* D_addr ((k * 128) + j) 0 */
int calAddrD_addr0( int i, int j, int k) {
    int result = (((k * 128) + j)) * 8 / 64;
    return result;
}
/* F_addr ((i * 128) + j) 1 */
int calAddrF_addr1( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* F_addr ((i * 128) + j) 2 */
int calAddrF_addr2( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* G_addr ((i * 128) + j) 0 */
int calAddrG_addr0( int i, int j) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* E_addr ((i * 128) + k) 3 */
int calAddrE_addr3( int i, int j, int k) {
    int result = (((i * 128) + k)) * 8 / 64;
    return result;
}
/* F_addr ((k * 128) + j) 3 */
int calAddrF_addr3( int i, int j, int k) {
    int result = (((k * 128) + j)) * 8 / 64;
    return result;
}
/* G_addr ((i * 128) + j) 1 */
int calAddrG_addr1( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* G_addr ((i * 128) + j) 2 */
int calAddrG_addr2( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
void ref_E_addr1() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[E_addr1]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[E_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrE_addr0( nv[nvi][0], nv[nvi][1]) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrE_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "E_addr1");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                    int kLB2 = 0;
                    if ( iLB0 == i_Start && jLB1 == j_Start ) {
                        kLB2 = k_Start;
                    }
                    for ( int k = kLB2; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[A_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[B_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr1");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 1 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr1");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 1 refNumber[LoopRefTree->AA]: 2 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = 0;
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            int iLB6 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB8 = 0;
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr1");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 1 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_E_addr2() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[E_addr2]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[E_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrE_addr0( nv[nvi][0], nv[nvi][1]) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrE_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "E_addr2");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 0 */
#endif
                    int kLB2 = 0;
                    if ( iLB0 == i_Start && jLB1 == j_Start ) {
                        kLB2 = k_Start;
                    }
                    for ( int k = kLB2; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[A_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[B_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr2");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr2");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 2 refNumber[LoopRefTree->AA]: 2 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = 0;
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            int iLB6 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB8 = 0;
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr2");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 2 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_E_addr0() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[E_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[E_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrE_addr0( nv[nvi][0], nv[nvi][1]) == calAddrE_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrE_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "E_addr0");
#endif
                                goto EndSample;
                            }
                        }
                        if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start) { cntStart = true; }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                    int kLB2 = 0;
                    for ( int k = kLB2; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[A_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[B_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 2 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = 0;
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            int iLB6 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB8 = 0;
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_C_addr0() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[C_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[F_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB5 = 0;
                    if ( iLB3 == i_Start && jLB4 == j_Start ) {
                        kLB5 = k_Start;
                    }
                    for ( int k = kLB5; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[C_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrC_addr0( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrC_addr0(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrC_addr0(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "C_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_D_addr0() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[D_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[F_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB5 = 0;
                    if ( iLB3 == i_Start && jLB4 == j_Start ) {
                        kLB5 = k_Start;
                    }
                    for ( int k = kLB5; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[C_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrD_addr0( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrD_addr0(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrD_addr0(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "D_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_A_addr0() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[A_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[E_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB2 = 0;
                    if ( iLB0 == i_Start && jLB1 == j_Start ) {
                        kLB2 = k_Start;
                    }
                    for ( int k = kLB2; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[A_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrA_addr0( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrA_addr0(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "A_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[B_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_B_addr0() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[B_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[E_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB2 = 0;
                    if ( iLB0 == i_Start && jLB1 == j_Start ) {
                        kLB2 = k_Start;
                    }
                    for ( int k = kLB2; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[A_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[B_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrB_addr0( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrB_addr0(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrB_addr0(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "B_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_G_addr0() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[G_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB6 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                if ( iLB6 == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrG_addr0( nv[nvi][0], nv[nvi][1]) == calAddrG_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrG_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "G_addr0");
#endif
                                goto EndSample;
                            }
                        }
                        if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start) { cntStart = true; }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                    int kLB8 = 0;
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrG_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrG_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrG_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "G_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrG_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrG_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrG_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "G_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 2 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_E_addr3() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[E_addr3]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB6 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                if ( iLB6 == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB8 = 0;
                    if ( iLB6 == i_Start && jLB7 == j_Start ) {
                        kLB8 = k_Start;
                    }
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrE_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrE_addr3(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrE_addr3(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "E_addr3");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 3 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_F_addr3() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[F_addr3]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB6 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                if ( iLB6 == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB8 = 0;
                    if ( iLB6 == i_Start && jLB7 == j_Start ) {
                        kLB8 = k_Start;
                    }
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr3(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr3(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr3");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 3 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_F_addr0() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[F_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[F_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrF_addr0( nv[nvi][0], nv[nvi][1]) == calAddrF_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrF_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "F_addr0");
#endif
                                goto EndSample;
                            }
                        }
                        if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start) { cntStart = true; }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                    int kLB5 = 0;
                    for ( int k = kLB5; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[C_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 2 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = 0;
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            int iLB6 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB8 = 0;
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr0");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 0 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_F_addr1() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[F_addr1]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[F_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrF_addr0( nv[nvi][0], nv[nvi][1]) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrF_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "F_addr1");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                    int kLB5 = 0;
                    if ( iLB3 == i_Start && jLB4 == j_Start ) {
                        kLB5 = k_Start;
                    }
                    for ( int k = kLB5; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[C_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr1");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 1 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr1");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 1 refNumber[LoopRefTree->AA]: 2 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = 0;
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            int iLB6 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB8 = 0;
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr1");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 1 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_F_addr2() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[F_addr2]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[F_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrF_addr0( nv[nvi][0], nv[nvi][1]) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrF_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "F_addr2");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 0 */
#endif
                    int kLB5 = 0;
                    if ( iLB3 == i_Start && jLB4 == j_Start ) {
                        kLB5 = k_Start;
                    }
                    for ( int k = kLB5; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[C_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr2");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr2");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 2 refNumber[LoopRefTree->AA]: 2 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = 0;
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            int iLB6 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB8 = 0;
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrF_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrF_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrF_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "F_addr2");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 2 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_G_addr1() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[G_addr1]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB6 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                if ( iLB6 == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrG_addr0( nv[nvi][0], nv[nvi][1]) == calAddrG_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrG_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "G_addr1");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                    int kLB8 = 0;
                    if ( iLB6 == i_Start && jLB7 == j_Start ) {
                        kLB8 = k_Start;
                    }
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrG_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrG_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrG_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "G_addr1");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 1 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrG_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrG_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrG_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "G_addr1");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 1 refNumber[LoopRefTree->AA]: 2 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_G_addr2() {
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2097;) {
SAMPLE:
        int i_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (128 - 0) + 0;
        if ( (128 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (128 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[G_addr2]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        auto BLIST = new int[THREAD_NUM][2];
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 3  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        int chunk_size, chunk_num, c_Start, ci_Start;
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 0) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB6 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB6 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB7 = 0;
                if ( iLB6 == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 128; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                    cout << "Generate interleaved iteration for (";
                    for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                        cout << *it;
                        if (it != v.end()) { cout << ", "; }
                    }
                    cout << ")" << endl;
#endif
                    for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                        vector<int> tmp;
                        for (int vi = 0; vi < v.size(); vi++ ) {
                            if (vi == 0) {
                                tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                            } else {
                                tmp.push_back(v[vi]);
                            }
                        }
                        if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                        cout << "(";
                        for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                            cout << *it << ", ";
                        }
                        cout << ")" << endl;
#endif
                    }
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[G_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrG_addr0( nv[nvi][0], nv[nvi][1]) == calAddrG_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrG_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "G_addr2");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 0 */
#endif
                    int kLB8 = 0;
                    if ( iLB6 == i_Start && jLB7 == j_Start ) {
                        kLB8 = k_Start;
                    }
                    for ( int k = kLB8; k < 128; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ", " << k << ", "")" << endl;
#endif
                        vector<int> v = { i, j, k };
                        /* Interleaving */
                        t_Start = ((i - 0) / chunk_size) % THREAD_NUM;
#ifdef DEBUG
                        cout << "Generate interleaved iteration for (";
                        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
                            cout << *it;
                            if (it != v.end()) { cout << ", "; }
                        }
                        cout << ")" << endl;
#endif
                        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {
                            vector<int> tmp;
                            for (int vi = 0; vi < v.size(); vi++ ) {
                                if (vi == 0) {
                                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));
                                } else {
                                    tmp.push_back(v[vi]);
                                }
                            }
                            if (tmp.size() > 0) { nv[tid] = tmp; }
#ifdef DEBUG
                            cout << "(";
                            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {
                                cout << *it << ", ";
                            }
                            cout << ")" << endl;
#endif
                        }
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[E_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[F_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                            }
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrG_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrG_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrG_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "G_addr2");
#endif
                                    goto EndSample;
                                }
                            }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[G_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrG_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrG_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrG_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(RT, cnt, 1.0);
#else
                                    refRTHistoCal(refRT, cnt, 1.0, "G_addr2");
#endif
                                    goto EndSample;
                                }
                            }
                            if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start && nv[nvi][2] == k_Start) { cntStart = true; }
                            }
#ifdef DEBUG
                            cout << endl;
                            /* useID: 2 refNumber[LoopRefTree->AA]: 2 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
int main() {
    /* 15 */
    std::thread t_E_addr_1(ref_E_addr1);
    std::thread t_E_addr_2(ref_E_addr2);
    std::thread t_E_addr_0(ref_E_addr0);
    std::thread t_C_addr_0(ref_C_addr0);
    std::thread t_D_addr_0(ref_D_addr0);
    std::thread t_A_addr_0(ref_A_addr0);
    std::thread t_B_addr_0(ref_B_addr0);
    std::thread t_G_addr_0(ref_G_addr0);
    std::thread t_E_addr_3(ref_E_addr3);
    std::thread t_F_addr_3(ref_F_addr3);
    std::thread t_F_addr_0(ref_F_addr0);
    std::thread t_F_addr_1(ref_F_addr1);
    std::thread t_F_addr_2(ref_F_addr2);
    std::thread t_G_addr_1(ref_G_addr1);
    std::thread t_G_addr_2(ref_G_addr2);
    t_E_addr_1.join();
    t_E_addr_2.join();
    t_E_addr_0.join();
    t_C_addr_0.join();
    t_D_addr_0.join();
    t_A_addr_0.join();
    t_B_addr_0.join();
    t_G_addr_0.join();
    t_E_addr_3.join();
    t_F_addr_3.join();
    t_F_addr_0.join();
    t_F_addr_1.join();
    t_F_addr_2.join();
    t_G_addr_1.join();
    t_G_addr_2.join();
    group_uniform_smoothing();
    rtMerge();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: mm2 */ 
