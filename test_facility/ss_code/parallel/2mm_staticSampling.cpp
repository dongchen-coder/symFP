
 /* Start to analysis array index
Array index info: Total number of references: 11
tmp.addr ((i * 4) + j)
A.addr ((i * 4) + k)
C.addr ((k * 4) + j)
D.addr ((i * 4) + j)
D.addr ((i * 4) + j)
B.addr ((k * 4) + j)
tmp.addr ((i * 4) + j)
tmp.addr ((i * 4) + j)
D.addr ((i * 4) + j)
D.addr ((i * 4) + j)
tmp.addr ((i * 4) + k)

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
--Loop Bound: (0, 4)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 4)
----Loop inc: (j + 1)
----Loop predicate: <
------array access tmp.addr ((i * 4) + j)
------k
------Loop Bound: (0, 4)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 4) + k)
--------array access B.addr ((k * 4) + j)
--------array access tmp.addr ((i * 4) + j)
--------array access tmp.addr ((i * 4) + j)
--i
--Loop Bound: (0, 4)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 4)
----Loop inc: (j + 1)
----Loop predicate: <
------array access D.addr ((i * 4) + j)
------array access D.addr ((i * 4) + j)
------k
------Loop Bound: (0, 4)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access tmp.addr ((i * 4) + k)
--------array access C.addr ((k * 4) + j)
--------array access D.addr ((i * 4) + j)
--------array access D.addr ((i * 4) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 4
------Sample number: 16
--------Sample number: 64
----Sample number: 4
------Sample number: 16
--------Sample number: 64
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
#include <map>
#include <set>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>
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
std::map<uint64_t, double> RT;
std::map<uint64_t, double> MR;
void rtHistoCal( int rt, int val ) {
    if ( val <= 0) {
;        return;
    }
    if (RT.find(rt) == RT.end()) { 
        RT[rt] = val;
    } else {
        RT[rt] += val;
    }
    return;
}
void subBlkRT(int rt) {
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
                rtHistoCal(b - diff, 1);
                break;
            }
        }
    }
    else {
        rtHistoCal(pow(2, msb-1), 1);
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
        cout << it->first << ", " << it->second << "\n";
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
/* tmp_addr ((i * 4) + j) 0 */
int calAddrtmp_addr0( int i, int j) {
    int result = (((i * 4) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 4) + k) 0 */
int calAddrA_addr0( int i, int j, int k) {
    int result = (((i * 4) + k)) * 8 / 64;
    return result;
}
/* B_addr ((k * 4) + j) 0 */
int calAddrB_addr0( int i, int j, int k) {
    int result = (((k * 4) + j)) * 8 / 64;
    return result;
}
/* tmp_addr ((i * 4) + j) 1 */
int calAddrtmp_addr1( int i, int j, int k) {
    int result = (((i * 4) + j)) * 8 / 64;
    return result;
}
/* tmp_addr ((i * 4) + j) 2 */
int calAddrtmp_addr2( int i, int j, int k) {
    int result = (((i * 4) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 4) + j) 0 */
int calAddrD_addr0( int i, int j) {
    int result = (((i * 4) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 4) + j) 1 */
int calAddrD_addr1( int i, int j) {
    int result = (((i * 4) + j)) * 8 / 64;
    return result;
}
/* tmp_addr ((i * 4) + k) 3 */
int calAddrtmp_addr3( int i, int j, int k) {
    int result = (((i * 4) + k)) * 8 / 64;
    return result;
}
/* C_addr ((k * 4) + j) 0 */
int calAddrC_addr0( int i, int j, int k) {
    int result = (((k * 4) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 4) + j) 2 */
int calAddrD_addr2( int i, int j, int k) {
    int result = (((i * 4) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 4) + j) 3 */
int calAddrD_addr3( int i, int j, int k) {
    int result = (((i * 4) + j)) * 8 / 64;
    return result;
}
void ref_tmp_addr0() {
/* for (i, 0, 4) */
/* for (j, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[tmp_addr0]Samples: " << idx_string << endl;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[tmp_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrtmp_addr0( nv[nvi][0], nv[nvi][1]) == calAddrtmp_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrtmp_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
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
                    for ( int k = kLB2; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[tmp_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            int iLB3 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                for ( int j = jLB4; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[D_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB5 = 0;
                    for ( int k = kLB5; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[D_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[D_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
/* for (i, 0, 4) */
/* for (j, 0, 4) */
/* for (k, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 64;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (4 - 0) + 0;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[tmp_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                    for ( int k = kLB2; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[tmp_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[tmp_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
/* for (i, 0, 4) */
/* for (j, 0, 4) */
/* for (k, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 64;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (4 - 0) + 0;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[D_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                    for ( int k = kLB5; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[C_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrC_addr0( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrC_addr0(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrC_addr0(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[D_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[D_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
void ref_D_addr2() {
/* for (i, 0, 4) */
/* for (j, 0, 4) */
/* for (k, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 64;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (4 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[D_addr2]Samples: " << idx_string << endl;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrD_addr0( nv[nvi][0], nv[nvi][1]) == calAddrD_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrD_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
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
                            cout  << "[D_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrD_addr1( nv[nvi][0], nv[nvi][1]) == calAddrD_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrD_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
#endif
                                goto EndSample;
                            }
                        }
                    }
#ifdef DEBUG
                    cout << endl;
                    /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
#endif
                    int kLB5 = 0;
                    if ( iLB3 == i_Start && jLB4 == j_Start ) {
                        kLB5 = k_Start;
                    }
                    for ( int k = kLB5; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[D_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrD_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrD_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrD_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                        /* iterate thread local iteration space mapping code after interleaving */
                        for (int nvi = 0; nvi < nv.size(); nvi++) {
                            /* Remove those invalid interleaving */
                            if (nv[nvi].size() <= 0) { continue; }
                            if (nv[nvi][0] > BLIST[nvi][1]) { break; }
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout  << "[D_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrD_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrD_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrD_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
#endif
                                    goto EndSample;
                                }
                            }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 3 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_D_addr3() {
/* for (i, 0, 4) */
/* for (j, 0, 4) */
/* for (k, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 64;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (4 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[D_addr3]Samples: " << idx_string << endl;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrD_addr0( nv[nvi][0], nv[nvi][1]) == calAddrD_addr3(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrD_addr3(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
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
                            cout  << "[D_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrD_addr1( nv[nvi][0], nv[nvi][1]) == calAddrD_addr3(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrD_addr3(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
#endif
                                goto EndSample;
                            }
                        }
                    }
#ifdef DEBUG
                    cout << endl;
                    /* useID: 3 refNumber[LoopRefTree->AA]: 1 */
#endif
                    int kLB5 = 0;
                    if ( iLB3 == i_Start && jLB4 == j_Start ) {
                        kLB5 = k_Start;
                    }
                    for ( int k = kLB5; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[D_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrD_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrD_addr3(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrD_addr3(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[D_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrD_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrD_addr3(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrD_addr3(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_B_addr0() {
/* for (i, 0, 4) */
/* for (j, 0, 4) */
/* for (k, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 64;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (4 - 0) + 0;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[tmp_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                    for ( int k = kLB2; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[tmp_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[tmp_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
void ref_tmp_addr1() {
/* for (i, 0, 4) */
/* for (j, 0, 4) */
/* for (k, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 64;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (4 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[tmp_addr1]Samples: " << idx_string << endl;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[tmp_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrtmp_addr0( nv[nvi][0], nv[nvi][1]) == calAddrtmp_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrtmp_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
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
                    for ( int k = kLB2; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[tmp_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            int iLB3 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                for ( int j = jLB4; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[D_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB5 = 0;
                    for ( int k = kLB5; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr1(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr1(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[D_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[D_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
void ref_tmp_addr2() {
/* for (i, 0, 4) */
/* for (j, 0, 4) */
/* for (k, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 64;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (4 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[tmp_addr2]Samples: " << idx_string << endl;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB0 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB1 = 0;
                if ( iLB0 == i_Start ) {
                    jLB1 = j_Start;
                }
                for ( int j = jLB1; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[tmp_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrtmp_addr0( nv[nvi][0], nv[nvi][1]) == calAddrtmp_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrtmp_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
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
                    for ( int k = kLB2; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr1( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[tmp_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            int iLB3 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                for ( int j = jLB4; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[D_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                        }
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
#endif
                    int kLB5 = 0;
                    for ( int k = kLB5; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr2(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr2(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[D_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[D_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
/* for (i, 0, 4) */
/* for (j, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrD_addr0( nv[nvi][0], nv[nvi][1]) == calAddrD_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrD_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
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
                            cout  << "[D_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrD_addr1( nv[nvi][0], nv[nvi][1]) == calAddrD_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrD_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
#endif
                                goto EndSample;
                            }
                        }
                    }
#ifdef DEBUG
                    cout << endl;
                    /* useID: 0 refNumber[LoopRefTree->AA]: 1 */
#endif
                    int kLB5 = 0;
                    for ( int k = kLB5; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[D_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrD_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrD_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrD_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[D_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrD_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrD_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrD_addr0(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
#endif
                                    goto EndSample;
                                }
                            }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 3 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_D_addr1() {
/* for (i, 0, 4) */
/* for (j, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[D_addr1]Samples: " << idx_string << endl;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrD_addr0( nv[nvi][0], nv[nvi][1]) == calAddrD_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrD_addr1(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
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
                            cout  << "[D_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
#endif
                            if ( calAddrD_addr1( nv[nvi][0], nv[nvi][1]) == calAddrD_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                                cout << "[REUSE FIND] @ (" << calAddrD_addr1(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << "), " << cnt << ") " << endl;
                                rtHistoCal(cnt, 1);
#else
                                subBlkRT(cnt);
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
                    int kLB5 = 0;
                    for ( int k = kLB5; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[D_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrD_addr2( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrD_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrD_addr1(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[D_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrD_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrD_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrD_addr1(nv[nvi][0], nv[nvi][1]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
#endif
                                    goto EndSample;
                                }
                            }
                        }
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 3 */
#endif
                } // end of inner for loops
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_tmp_addr3() {
/* for (i, 0, 4) */
/* for (j, 0, 4) */
/* for (k, 0, 4) */
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 64;) {
SAMPLE:
        int i_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4 - 0) + 0;
        if ( (4 - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (4 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[tmp_addr3]Samples: " << idx_string << endl;
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
        chunk_num = (4 - 0) % (THREAD_NUM * chunk_size) == 0 ? (4 - 0) / (THREAD_NUM * chunk_size) : (4 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (4 - 0) / THREAD_NUM;
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
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            int iLB3 = 0;
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB3 = cid * (THREAD_NUM * chunk_size) + ci;
                }
                int jLB4 = 0;
                if ( iLB3 == i_Start ) {
                    jLB4 = j_Start;
                }
                for ( int j = jLB4; j < 4; j++) {
                    int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                    if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                    cout << "Iterate (" << i << ")" << endl;
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
                            cout  << "[D_addr0]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                            cout  << "[D_addr1]" << nv[nvi][0] << ", " << nv[nvi][1] << ", cnt: " << cnt << ")	";
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
                    for ( int k = kLB5; k < 4; k++) {
                        int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
                        if(i > BLIST[0][1]) { goto EndSample; }
#ifdef DEBUG
                        cout << "Iterate (" << i << ", " << j << ")" << endl;
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
                                cout  << "[tmp_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
#endif
                                if ( calAddrtmp_addr3( nv[nvi][0], nv[nvi][1], nv[nvi][2]) == calAddrtmp_addr3(i_Start, j_Start, k_Start)) {
#ifdef DEBUG
                                    cout << "[REUSE FIND] @ (" << calAddrtmp_addr3(nv[nvi][0], nv[nvi][1], nv[nvi][2]) << ", " << "(" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << "), " << cnt << ") " << endl;
                                    rtHistoCal(cnt, 1);
#else
                                    subBlkRT(cnt);
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
                                cout  << "[D_addr2]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
                                cout  << "[D_addr3]" << nv[nvi][0] << ", " << nv[nvi][1] << ", " << nv[nvi][2] << ", cnt: " << cnt << ")	";
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
int main() {
    ref_tmp_addr0();
    ref_A_addr0();
    ref_C_addr0();
    ref_D_addr2();
    ref_D_addr3();
    ref_B_addr0();
    ref_tmp_addr1();
    ref_tmp_addr2();
    ref_D_addr0();
    ref_D_addr1();
    ref_tmp_addr3();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: mm2 */ 
