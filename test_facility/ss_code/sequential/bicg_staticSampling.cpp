
 /* Start to analysis array index
Array index info: Total number of references: 9
q.addr i
A.addr ((i * 4096) + j)
p.addr j
q.addr i
s.addr j
r.addr i
A.addr ((i * 4096) + j)
s.addr j
q.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %nx
i32 %ny
double* %A
double* %r
double* %s
double* %p
double* %q

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 4096)
--Loop inc: (i + 1)
--Loop predicate: <
----array access q.addr i
----j
----Loop Bound: (0, 4096)
----Loop inc: (j + 1)
----Loop predicate: <
------array access s.addr j
------array access r.addr i
------array access A.addr ((i * 4096) + j)
------array access s.addr j
------array access q.addr i
------array access A.addr ((i * 4096) + j)
------array access p.addr j
------array access q.addr i

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 4096
------Sample number: 16777216
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
/* q_addr i 0 */
int calAddrq_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* s_addr j 0 */
int calAddrs_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* r_addr i 0 */
int calAddrr_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* A_addr ((i * 4096) + j) 0 */
int calAddrA_addr0( int i, int j) {
    int result = (((i * 4096) + j)) * 8 / 64;
    return result;
}
/* s_addr j 1 */
int calAddrs_addr1( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* q_addr i 1 */
int calAddrq_addr1( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* A_addr ((i * 4096) + j) 1 */
int calAddrA_addr1( int i, int j) {
    int result = (((i * 4096) + j)) * 8 / 64;
    return result;
}
/* p_addr j 0 */
int calAddrp_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* q_addr i 2 */
int calAddrq_addr2( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
void ref_q_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 4096;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[q_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 1  */
        /* Sampled IDV: i  */
        int iLB0 = i_Start;
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        int chunk_size = CHUNK_SIZE;
        int chunk_num = (4096 - 0) / (THREAD_NUM * chunk_size);
#else
        int chunk_size = (4096 - 0) / THREAD_NUM;
        int chunk_num = 1;
#endif
        /* Compute the number of chunks */
        int c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4096) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            int ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
#ifdef DEBUG
            cout << "Iterate (" << i << ")" << endl;
#endif
            vector<int> v = { i};
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
            /* Generating thread local iteration space mapping code */
            /* iterate thread local iteration space mapping code after interleaving */
            for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                /* Remove those invalid interleaving */
                if (it->size() <= 0) { continue; }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout  << "[q_addr0]" << it->at(0) << ", cnt: " << cnt << ")	";
#endif
                    if ( calAddrq_addr0( it->at(0)) == calAddrq_addr0(i_Start)) {
#ifdef DEBUG
                        cout << "[REUSE FIND] @ (" << calAddrq_addr0(it->at(0)) << ", " << "(" << it->at(0) << "), " << cnt << ") " << endl;
                        rtHistoCal(cnt, 1);
#else
                        subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
                if (it->at(0) == i_Start            ) { cntStart = true; }
            }
#ifdef DEBUG
            cout << endl;
            /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
            /* Generating thread local iteration space mapping code */
            for ( int j = jLB1; j < 4096; j++) {
                int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
#ifdef DEBUG
                cout << "Iterate (" << i << ")" << endl;
#endif
                vector<int> v = { i};
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[r_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrq_addr1( it->at(0), it->at(1)) == calAddrq_addr0(i_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrq_addr0(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[p_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr2]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrq_addr2( it->at(0), it->at(1)) == calAddrq_addr0(i_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrq_addr0(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
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
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_A_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16777216;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[A_addr1]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        int iLB0 = i_Start;
        int jLB1 = 0;
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        int chunk_size = CHUNK_SIZE;
        int chunk_num = (4096 - 0) / (THREAD_NUM * chunk_size);
#else
        int chunk_size = (4096 - 0) / THREAD_NUM;
        int chunk_num = 1;
#endif
        /* Compute the number of chunks */
        int c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4096) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            int ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
#ifdef DEBUG
            cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
            vector<int> v = { i, j};
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
            /* Generating thread local iteration space mapping code */
            /* iterate thread local iteration space mapping code after interleaving */
            for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                if (it->size() <= 0) { continue; }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout  << "[q_addr0]" << it->at(0) << ", cnt: " << cnt << ")	";
#endif
                }
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
#endif
            /* Generating thread local iteration space mapping code */
                if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j++) {
                int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
#ifdef DEBUG
                cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
                vector<int> v = { i, j};
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[r_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrA_addr0( it->at(0), it->at(1)) == calAddrA_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrA_addr1(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrA_addr1( it->at(0), it->at(1)) == calAddrA_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrA_addr1(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                            goto EndSample;
                        }
                    }
                    if (it->at(0) == i_Start && it->at(1) == j_Start                ) { cntStart = true; }
                }
#ifdef DEBUG
                cout << endl;
                /* useID: 1 refNumber[LoopRefTree->AA]: 1 */
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[p_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr2]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
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
void ref_p_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16777216;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[p_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        int iLB0 = i_Start;
        int jLB1 = 0;
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        int chunk_size = CHUNK_SIZE;
        int chunk_num = (4096 - 0) / (THREAD_NUM * chunk_size);
#else
        int chunk_size = (4096 - 0) / THREAD_NUM;
        int chunk_num = 1;
#endif
        /* Compute the number of chunks */
        int c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4096) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            int ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
#ifdef DEBUG
            cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
            vector<int> v = { i, j};
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
            /* Generating thread local iteration space mapping code */
            /* iterate thread local iteration space mapping code after interleaving */
            for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                if (it->size() <= 0) { continue; }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout  << "[q_addr0]" << it->at(0) << ", cnt: " << cnt << ")	";
#endif
                }
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
#endif
            /* Generating thread local iteration space mapping code */
                if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j++) {
                int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
#ifdef DEBUG
                cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
                vector<int> v = { i, j};
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[r_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[p_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrp_addr0( it->at(0), it->at(1)) == calAddrp_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrp_addr0(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                            goto EndSample;
                        }
                    }
                    if (it->at(0) == i_Start && it->at(1) == j_Start                ) { cntStart = true; }
                }
#ifdef DEBUG
                cout << endl;
                /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr2]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
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
void ref_q_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16777216;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[q_addr2]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        int iLB0 = i_Start;
        int jLB1 = 0;
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        int chunk_size = CHUNK_SIZE;
        int chunk_num = (4096 - 0) / (THREAD_NUM * chunk_size);
#else
        int chunk_size = (4096 - 0) / THREAD_NUM;
        int chunk_num = 1;
#endif
        /* Compute the number of chunks */
        int c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4096) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            int ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
#ifdef DEBUG
            cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
            vector<int> v = { i, j};
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
            /* Generating thread local iteration space mapping code */
            /* iterate thread local iteration space mapping code after interleaving */
            for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                /* Remove those invalid interleaving */
                if (it->size() <= 0) { continue; }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout  << "[q_addr0]" << it->at(0) << ", cnt: " << cnt << ")	";
#endif
                    if ( calAddrq_addr0( it->at(0)) == calAddrq_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                        cout << "[REUSE FIND] @ (" << calAddrq_addr2(it->at(0)) << ", " << "(" << it->at(0) << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
                if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j++) {
                int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
#ifdef DEBUG
                cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
                vector<int> v = { i, j};
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[r_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrq_addr1( it->at(0), it->at(1)) == calAddrq_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrq_addr2(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[p_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr2]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrq_addr2( it->at(0), it->at(1)) == calAddrq_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrq_addr2(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                            goto EndSample;
                        }
                    }
                    if (it->at(0) == i_Start && it->at(1) == j_Start                ) { cntStart = true; }
                }
#ifdef DEBUG
                cout << endl;
                /* useID: 2 refNumber[LoopRefTree->AA]: 2 */
#endif
            } // end of inner for loops
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
void ref_s_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16777216;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[s_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        int iLB0 = i_Start;
        int jLB1 = 0;
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        int chunk_size = CHUNK_SIZE;
        int chunk_num = (4096 - 0) / (THREAD_NUM * chunk_size);
#else
        int chunk_size = (4096 - 0) / THREAD_NUM;
        int chunk_num = 1;
#endif
        /* Compute the number of chunks */
        int c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4096) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            int ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
#ifdef DEBUG
            cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
            vector<int> v = { i, j};
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
            /* Generating thread local iteration space mapping code */
            /* iterate thread local iteration space mapping code after interleaving */
            for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                if (it->size() <= 0) { continue; }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout  << "[q_addr0]" << it->at(0) << ", cnt: " << cnt << ")	";
#endif
                }
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
#endif
            /* Generating thread local iteration space mapping code */
                if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j++) {
                int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
#ifdef DEBUG
                cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
                vector<int> v = { i, j};
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrs_addr0( it->at(0), it->at(1)) == calAddrs_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrs_addr0(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                            goto EndSample;
                        }
                    }
                    if (it->at(0) == i_Start && it->at(1) == j_Start                ) { cntStart = true; }
                }
#ifdef DEBUG
                cout << endl;
                /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[r_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrs_addr1( it->at(0), it->at(1)) == calAddrs_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrs_addr0(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[p_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr2]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
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
void ref_r_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16777216;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[r_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        int iLB0 = i_Start;
        int jLB1 = 0;
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        int chunk_size = CHUNK_SIZE;
        int chunk_num = (4096 - 0) / (THREAD_NUM * chunk_size);
#else
        int chunk_size = (4096 - 0) / THREAD_NUM;
        int chunk_num = 1;
#endif
        /* Compute the number of chunks */
        int c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4096) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            int ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
#ifdef DEBUG
            cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
            vector<int> v = { i, j};
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
            /* Generating thread local iteration space mapping code */
            /* iterate thread local iteration space mapping code after interleaving */
            for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                if (it->size() <= 0) { continue; }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout  << "[q_addr0]" << it->at(0) << ", cnt: " << cnt << ")	";
#endif
                }
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
#endif
            /* Generating thread local iteration space mapping code */
                if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j++) {
                int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
#ifdef DEBUG
                cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
                vector<int> v = { i, j};
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[r_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrr_addr0( it->at(0), it->at(1)) == calAddrr_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrr_addr0(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                            goto EndSample;
                        }
                    }
                    if (it->at(0) == i_Start && it->at(1) == j_Start                ) { cntStart = true; }
                }
#ifdef DEBUG
                cout << endl;
                /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[p_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr2]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
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
void ref_A_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16777216;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[A_addr0]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        int iLB0 = i_Start;
        int jLB1 = 0;
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        int chunk_size = CHUNK_SIZE;
        int chunk_num = (4096 - 0) / (THREAD_NUM * chunk_size);
#else
        int chunk_size = (4096 - 0) / THREAD_NUM;
        int chunk_num = 1;
#endif
        /* Compute the number of chunks */
        int c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4096) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            int ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
#ifdef DEBUG
            cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
            vector<int> v = { i, j};
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
            /* Generating thread local iteration space mapping code */
            /* iterate thread local iteration space mapping code after interleaving */
            for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                if (it->size() <= 0) { continue; }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout  << "[q_addr0]" << it->at(0) << ", cnt: " << cnt << ")	";
#endif
                }
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
#endif
            /* Generating thread local iteration space mapping code */
                if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j++) {
                int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
#ifdef DEBUG
                cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
                vector<int> v = { i, j};
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[r_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrA_addr0( it->at(0), it->at(1)) == calAddrA_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrA_addr0(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                            goto EndSample;
                        }
                    }
                    if (it->at(0) == i_Start && it->at(1) == j_Start                ) { cntStart = true; }
                }
#ifdef DEBUG
                cout << endl;
                /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrA_addr1( it->at(0), it->at(1)) == calAddrA_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrA_addr0(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[p_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr2]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
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
void ref_s_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16777216;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[s_addr1]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        int iLB0 = i_Start;
        int jLB1 = 0;
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        int chunk_size = CHUNK_SIZE;
        int chunk_num = (4096 - 0) / (THREAD_NUM * chunk_size);
#else
        int chunk_size = (4096 - 0) / THREAD_NUM;
        int chunk_num = 1;
#endif
        /* Compute the number of chunks */
        int c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4096) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            int ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
#ifdef DEBUG
            cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
            vector<int> v = { i, j};
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
            /* Generating thread local iteration space mapping code */
            /* iterate thread local iteration space mapping code after interleaving */
            for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                if (it->size() <= 0) { continue; }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout  << "[q_addr0]" << it->at(0) << ", cnt: " << cnt << ")	";
#endif
                }
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
#endif
            /* Generating thread local iteration space mapping code */
                if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j++) {
                int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
#ifdef DEBUG
                cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
                vector<int> v = { i, j};
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrs_addr0( it->at(0), it->at(1)) == calAddrs_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrs_addr1(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[r_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrs_addr1( it->at(0), it->at(1)) == calAddrs_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrs_addr1(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                            goto EndSample;
                        }
                    }
                    if (it->at(0) == i_Start && it->at(1) == j_Start                ) { cntStart = true; }
                }
#ifdef DEBUG
                cout << endl;
                /* useID: 1 refNumber[LoopRefTree->AA]: 1 */
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[p_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr2]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
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
void ref_q_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16777216;) {
SAMPLE:
        int i_Start = rand() % (4096 - 0) + 0;
        if ( (4096 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (4096 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
#ifdef DEBUG
        cout << "[q_addr1]Samples: " << idx_string << endl;
#endif
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int t_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 2  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        int iLB0 = i_Start;
        int jLB1 = 0;
        /* Vector that contains the interleaved iteration, avoid duplicate declaration */
        vector<vector<int>> nv(THREAD_NUM);
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        int chunk_size = CHUNK_SIZE;
        int chunk_num = (4096 - 0) / (THREAD_NUM * chunk_size);
#else
        int chunk_size = (4096 - 0) / THREAD_NUM;
        int chunk_num = 1;
#endif
        /* Compute the number of chunks */
        int c_Start = (i_Start - 1) / (THREAD_NUM * chunk_size);
        cout << "c_Start = " << c_Start << ", chunk_num = " << chunk_num << endl;
        for (int cid = c_Start; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 4096) - 1;
#ifdef DEBUG
                cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            /* Iterate within a chunk */
            int ci_Start = 0;
            if (cid == c_Start) {
                ci_Start = (i_Start - 0) % chunk_size;
            }
            for ( int ci = ci_Start; ci < chunk_size; ci++) {
                if ( cid != c_Start || ci != ci_Start ) {
                    iLB0 = cid * (THREAD_NUM * chunk_size) + ci;
                }
#ifdef DEBUG
            cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
            vector<int> v = { i, j};
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
            /* Generating thread local iteration space mapping code */
            /* iterate thread local iteration space mapping code after interleaving */
            for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                /* Remove those invalid interleaving */
                if (it->size() <= 0) { continue; }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout  << "[q_addr0]" << it->at(0) << ", cnt: " << cnt << ")	";
#endif
                    if ( calAddrq_addr0( it->at(0)) == calAddrq_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                        cout << "[REUSE FIND] @ (" << calAddrq_addr1(it->at(0)) << ", " << "(" << it->at(0) << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
                if ( iLB0 == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 4096; j++) {
                int i = cid * (THREAD_NUM * chunk_size) + ci + 0;
#ifdef DEBUG
                cout << "Iterate (" << i << ", " << j << ")" << endl;
#endif
                vector<int> v = { i, j};
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
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[r_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[s_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrq_addr1( it->at(0), it->at(1)) == calAddrq_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrq_addr1(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                            goto EndSample;
                        }
                    }
                    if (it->at(0) == i_Start && it->at(1) == j_Start                ) { cntStart = true; }
                }
#ifdef DEBUG
                cout << endl;
                /* useID: 1 refNumber[LoopRefTree->AA]: 1 */
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[A_addr1]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[p_addr0]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                    }
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
#endif
                /* Generating thread local iteration space mapping code */
                /* iterate thread local iteration space mapping code after interleaving */
                for (vector<vector<int>>::iterator it = nv.begin(); it != nv.end(); it++) {
                    /* Remove those invalid interleaving */
                    if (it->size() <= 0) { continue; }
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout  << "[q_addr2]" << it->at(0) << ", " << it->at(1) << ", cnt: " << cnt << ")	";
#endif
                        if ( calAddrq_addr2( it->at(0), it->at(1)) == calAddrq_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            cout << "[REUSE FIND] @ (" << calAddrq_addr1(it->at(0), it->at(1)) << ", " << "(" << it->at(0) << ", " << it->at(1) << "), " << cnt << ") " << endl;
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
            } // end of outer for - ci loops
        } // end of outer for - cid loops
EndSample:
        s++;
        }
}
int main() {
    ref_q_addr0();
    ref_A_addr1();
    ref_p_addr0();
    ref_q_addr2();
    ref_s_addr0();
    ref_r_addr0();
    ref_A_addr0();
    ref_s_addr1();
    ref_q_addr1();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: bicg_cpu */ 
