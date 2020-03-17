
 /* Start to analysis array index
Array index info: Total number of references: 6
a.addr ((1026 * i) + j)
a.addr (((1026 * i) + j) + 1)
b.addr ((1026 * i) + j)
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
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (1, 1025)
----Loop inc: (j + 1)
----Loop predicate: <
------array access a.addr ((1026 * i) + j)
------array access a.addr (((1026 * i) + j) + 1)
------array access a.addr (((1026 * i) + j) - 1)
------array access a.addr ((1026 * (i - 1)) + j)
------array access a.addr ((1026 * (i + 1)) + j)
------array access b.addr ((1026 * i) + j)

Finish analysis loops */ 
/* # of Out-most Loops: 1 */ 

 /* Start transform loop tree
	for.cond
	for.cond1
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
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------
--------------------
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------
--------------------
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------
--------------------
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------
--------------------
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------

Finish transform loop tree */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 102
------Sample number: 10485
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
std::map<string, map<uint64_t, double>> refRT;
std::map<uint64_t, double> RT;
std::map<uint64_t, double> MR;
void rtHistoCal( map<uint64_t, double> &rth, uint64_t rt, double val ) {
    if ( val <= 0) {
;        return;
    }
    if (rth.find(rt) == rth.end()) { 
        rth[rt] = val;
    } else {
        rth[rt] += val;
    }
    return;
}

void refRTHistoCal(map<string, map<uint64_t, double>> &rth, uint64_t rt, double val, string ref ) {
    if ( val <= 0) {
;        return;
    }
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
/* Smoothing the per-reference RTHisto Uniformly. Equally splite the RT to a range of bins */
void uniform_smoothing(map<uint64_t, double> &rth, bool enable, bool scale) {
    map<uint64_t, double> tmp;
    double sum = accumulate(begin(rth), end(rth), 0.0, [](const double previous, const pair<uint64_t, double>& p) { return previous + p.second; });
    double sum_P = 0.0;
    for(map<uint64_t, double>::iterator it = rth.begin(); it != rth.end(); ++it) {
        uint64_t mu = it->first;
        sum_P += it->second;
        /* Do uniform distribution for all ri, distribute from ri / THREAD_NUM to ri * THREAD_NUM  */
        if (it->second <= 0.0) { continue; }
        if (mu < THREAD_NUM) { mu = mu * THREAD_NUM; }
        uint64_t start_b = (mu / THREAD_NUM) >= 1 ? mu / THREAD_NUM : 1;
        uint64_t end_b = mu * THREAD_NUM;
        double split_val = it->second / (end_b - start_b + 1);
        for (int b = start_b; b <= end_b; b++) {
            rtHistoCal(tmp, b, split_val);
        }
    }
    rth.clear();
    for(map<uint64_t, double>::iterator it = tmp.begin(); it != tmp.end(); ++it) {
        subBlkRT(rth, it->first, it->second);
    }
    return;
}
/* Array a_addr	i j */ 
/* Array a_addr	i j */ 
/* Array b_addr	i j */ 
/* Array a_addr	i j */ 
/* Array a_addr	i j */ 
/* Array a_addr	i j */ 
void group_uniform_smoothing() {
    for (map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {
    /* Cold Program - No scaling */
    uniform_smoothing(it->second, false, false);
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
/* a_addr ((1026 * i) + j) 0 */
int calAddra_addr0( int i, int j) {
    int result = (((1026 * i) + j)) * 8 / 64;
    return result;
}
/* a_addr (((1026 * i) + j) + 1) 1 */
int calAddra_addr1( int i, int j) {
    int result = ((((1026 * i) + j) + 1)) * 8 / 64;
    return result;
}
/* a_addr (((1026 * i) + j) - 1) 2 */
int calAddra_addr2( int i, int j) {
    int result = ((((1026 * i) + j) - 1)) * 8 / 64;
    return result;
}
/* a_addr ((1026 * (i - 1)) + j) 3 */
int calAddra_addr3( int i, int j) {
    int result = (((1026 * (i - 1)) + j)) * 8 / 64;
    return result;
}
/* a_addr ((1026 * (i + 1)) + j) 4 */
int calAddra_addr4( int i, int j) {
    int result = (((1026 * (i + 1)) + j)) * 8 / 64;
    return result;
}
/* b_addr ((1026 * i) + j) 0 */
int calAddrb_addr0( int i, int j) {
    int result = (((1026 * i) + j)) * 8 / 64;
    return result;
}
void ref_a_addr0() {
/* for (i, 1, 1025) */
/* for (j, 1, 1025) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10485;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[a_addr0]Samples: " << idx_string << endl;
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
        chunk_num = (1025 - 1) % (THREAD_NUM * chunk_size) == 0 ? (1025 - 1) / (THREAD_NUM * chunk_size) : (1025 - 1) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1025 - 1) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  1+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(1 + (cid * THREAD_NUM + t + 1) * chunk_size, 1025) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 1) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 1;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 1025; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 1;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 1) / chunk_size) % THREAD_NUM;
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
                            cout  << "[a_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr0( nv[nvi][0], nv[nvi][1]) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr0");
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
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr1( nv[nvi][0], nv[nvi][1]) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr0");
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
                            cout  << "[a_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr2( nv[nvi][0], nv[nvi][1]) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr0");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 2 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr3( nv[nvi][0], nv[nvi][1]) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr0");
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
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr4]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr4( nv[nvi][0], nv[nvi][1]) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr0");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 4 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[b_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_a_addr1() {
/* for (i, 1, 1025) */
/* for (j, 1, 1025) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10485;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[a_addr1]Samples: " << idx_string << endl;
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
        chunk_num = (1025 - 1) % (THREAD_NUM * chunk_size) == 0 ? (1025 - 1) / (THREAD_NUM * chunk_size) : (1025 - 1) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1025 - 1) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  1+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(1 + (cid * THREAD_NUM + t + 1) * chunk_size, 1025) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 1) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 1;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 1025; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 1;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 1) / chunk_size) % THREAD_NUM;
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
                            cout  << "[a_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr0( nv[nvi][0], nv[nvi][1]) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr1(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr1");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr1( nv[nvi][0], nv[nvi][1]) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr1(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr1");
#endif
                                goto EndSample;
                            }
                        }
                        if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start) { cntStart = true; }
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
                            cout  << "[a_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr2( nv[nvi][0], nv[nvi][1]) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr1(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr1");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 2 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr3( nv[nvi][0], nv[nvi][1]) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr1(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr1");
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
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr4]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr4( nv[nvi][0], nv[nvi][1]) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr1(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr1");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 4 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[b_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_b_addr0() {
/* for (i, 1, 1025) */
/* for (j, 1, 1025) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10485;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[b_addr0]Samples: " << idx_string << endl;
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
        chunk_num = (1025 - 1) % (THREAD_NUM * chunk_size) == 0 ? (1025 - 1) / (THREAD_NUM * chunk_size) : (1025 - 1) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1025 - 1) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  1+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(1 + (cid * THREAD_NUM + t + 1) * chunk_size, 1025) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 1) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 1;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 1025; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 1;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 1) / chunk_size) % THREAD_NUM;
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
                            cout  << "[a_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[a_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[a_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[a_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[a_addr4]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[b_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrb_addr0( nv[nvi][0], nv[nvi][1]) == calAddrb_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrb_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "b_addr0");
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
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_a_addr2() {
/* for (i, 1, 1025) */
/* for (j, 1, 1025) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10485;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[a_addr2]Samples: " << idx_string << endl;
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
        chunk_num = (1025 - 1) % (THREAD_NUM * chunk_size) == 0 ? (1025 - 1) / (THREAD_NUM * chunk_size) : (1025 - 1) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1025 - 1) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  1+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(1 + (cid * THREAD_NUM + t + 1) * chunk_size, 1025) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 1) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 1;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 1025; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 1;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 1) / chunk_size) % THREAD_NUM;
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
                            cout  << "[a_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr0( nv[nvi][0], nv[nvi][1]) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr2(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr2");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 0 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr1( nv[nvi][0], nv[nvi][1]) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr2(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr2");
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
                            cout  << "[a_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr2( nv[nvi][0], nv[nvi][1]) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr2(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr2");
#endif
                                goto EndSample;
                            }
                        }
                        if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start) { cntStart = true; }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 2 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr3( nv[nvi][0], nv[nvi][1]) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr2(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr2");
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
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr4]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr4( nv[nvi][0], nv[nvi][1]) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr2(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr2");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 4 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[b_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_a_addr3() {
/* for (i, 1, 1025) */
/* for (j, 1, 1025) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10485;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[a_addr3]Samples: " << idx_string << endl;
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
        chunk_num = (1025 - 1) % (THREAD_NUM * chunk_size) == 0 ? (1025 - 1) / (THREAD_NUM * chunk_size) : (1025 - 1) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1025 - 1) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  1+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(1 + (cid * THREAD_NUM + t + 1) * chunk_size, 1025) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 1) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 1;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 1025; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 1;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 1) / chunk_size) % THREAD_NUM;
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
                            cout  << "[a_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr0( nv[nvi][0], nv[nvi][1]) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr3(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr3");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 0 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr1( nv[nvi][0], nv[nvi][1]) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr3(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr3");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 1 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr2( nv[nvi][0], nv[nvi][1]) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr3(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr3");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 2 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr3( nv[nvi][0], nv[nvi][1]) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr3(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr3");
#endif
                                goto EndSample;
                            }
                        }
                        if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start) { cntStart = true; }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 3 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr4]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr4( nv[nvi][0], nv[nvi][1]) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr3(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr3");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 4 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[b_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_a_addr4() {
/* for (i, 1, 1025) */
/* for (j, 1, 1025) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10485;) {
SAMPLE:
        int i_Start = rand() % (1025 - 1) + 1;
        if ( (1025 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[a_addr4]Samples: " << idx_string << endl;
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
        chunk_num = (1025 - 1) % (THREAD_NUM * chunk_size) == 0 ? (1025 - 1) / (THREAD_NUM * chunk_size) : (1025 - 1) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1025 - 1) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
#ifdef DEBUG
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
#endif
        /* Generating thread local iteration space mapping code */
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  1+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(1 + (cid * THREAD_NUM + t + 1) * chunk_size, 1025) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 1) % chunk_size;
            }
            int iLB0 = i_Start;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 1;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 1025; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 1;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ", " << j << ", "")" << endl;
#endif
                    vector<int> v = { i, j };
                    /* Interleaving */
                    t_Start = ((i - 1) / chunk_size) % THREAD_NUM;
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
                            cout  << "[a_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr0( nv[nvi][0], nv[nvi][1]) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr4(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr4");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 0 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr1( nv[nvi][0], nv[nvi][1]) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr4(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr4");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 1 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr2( nv[nvi][0], nv[nvi][1]) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr4(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr4");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 2 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr3( nv[nvi][0], nv[nvi][1]) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr4(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr4");
#endif
                                goto EndSample;
                            }
                        }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 3 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        /* Remove those invalid interleaving */
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[a_addr4]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddra_addr4( nv[nvi][0], nv[nvi][1]) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddra_addr4(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(RT, cnt, 1.0);
#else
                                refRTHistoCal(refRT, cnt, 1.0, "a_addr4");
#endif
                                goto EndSample;
                            }
                        }
                        if (nv[nvi][0] == i_Start && nv[nvi][1] == j_Start) { cntStart = true; }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 4 */
#endif
                    /* iterate thread local iteration space mapping code after interleaving */
                    for (int nvi = 0; nvi < nv.size(); nvi++) {
                        if (nv[nvi].size() <= 0) { continue; }
                        if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout  << "[b_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
int main() {
    /* 6 */
    ref_a_addr0();
    ref_a_addr1();
    ref_b_addr0();
    ref_a_addr2();
    ref_a_addr3();
    ref_a_addr4();
    group_uniform_smoothing();
    rtMerge();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: stencil */ 
