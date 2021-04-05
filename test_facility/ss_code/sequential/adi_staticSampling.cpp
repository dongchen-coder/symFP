
 /* Start to analysis array index
Array index info: Total number of references: 34
v.addr ((0 + i) + 1)
p.addr (((i + 1) * 8192) + j)
u.addr (((i + 1) * 8192) + 0)
q.addr (((i + 1) * 8192) + 0)
u.addr (((j + 1) * 8192) + i)
p.addr (((i + 1) * 8192) + 0)
v.addr (0 + (i + 1))
q.addr (((i + 1) * 8192) + 0)
q.addr ((((i + 1) * 8192) + j) + 1)
v.addr (67100672 + (i + 1))
p.addr (((i + 1) * 8192) + (8190 - j))
p.addr ((((i + 1) * 8192) + j) + 1)
v.addr (((8190 - j) * 8192) + (i + 1))
u.addr (((i + 1) * 8192) + 0)
p.addr (((i + 1) * 8192) + 0)
u.addr ((((j + 1) * 8192) + i) + 1)
u.addr ((((j + 1) * 8192) + i) + 2)
q.addr (((i + 1) * 8192) + j)
p.addr (((i + 1) * 8192) + j)
q.addr ((((i + 1) * 8192) + j) + 1)
u.addr ((((i + 1) * 8192) + 8192) - 1)
p.addr (((((i + 1) * 8192) + 8192) - 2) - j)
u.addr (((((i + 1) * 8192) + 8192) - 1) - j)
q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
v.addr ((((8191 - j) * 8192) + i) + 1)
q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
p.addr (((i + 1) * 8192) + j)
p.addr ((((i + 1) * 8192) + j) + 1)
v.addr (((i * 8192) + j) + 1)
v.addr ((((i + 1) * 8192) + j) + 1)
v.addr ((((i + 2) * 8192) + j) + 1)
q.addr (((i + 1) * 8192) + j)
p.addr (((i + 1) * 8192) + j)
u.addr (((((i + 1) * 8192) + 8192) - 2) - j)
BC Array cost info: Total number of arrays: 4
p.addr 24
q.addr 20
u.addr 19
v.addr 19
BC Average cost: 2.050000e+01

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
--i
--Loop Bound: (0, 8190)
--Loop inc: (i + 1)
--Loop predicate: <
----array access v.addr ((0 + i) + 1)
----array access p.addr (((i + 1) * 8192) + 0)
----array access v.addr (0 + (i + 1))
----array access q.addr (((i + 1) * 8192) + 0)
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p.addr (((i + 1) * 8192) + j)
------array access p.addr ((((i + 1) * 8192) + j) + 1)
------array access u.addr (((j + 1) * 8192) + i)
------array access u.addr ((((j + 1) * 8192) + i) + 1)
------array access u.addr ((((j + 1) * 8192) + i) + 2)
------array access q.addr (((i + 1) * 8192) + j)
------array access p.addr (((i + 1) * 8192) + j)
------array access q.addr ((((i + 1) * 8192) + j) + 1)
----array access v.addr (67100672 + (i + 1))
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p.addr (((i + 1) * 8192) + (8190 - j))
------array access v.addr ((((8191 - j) * 8192) + i) + 1)
------array access q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
------array access v.addr (((8190 - j) * 8192) + (i + 1))
--i
--Loop Bound: (0, 8190)
--Loop inc: (i + 1)
--Loop predicate: <
----array access u.addr (((i + 1) * 8192) + 0)
----array access p.addr (((i + 1) * 8192) + 0)
----array access u.addr (((i + 1) * 8192) + 0)
----array access q.addr (((i + 1) * 8192) + 0)
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p.addr (((i + 1) * 8192) + j)
------array access p.addr ((((i + 1) * 8192) + j) + 1)
------array access v.addr (((i * 8192) + j) + 1)
------array access v.addr ((((i + 1) * 8192) + j) + 1)
------array access v.addr ((((i + 2) * 8192) + j) + 1)
------array access q.addr (((i + 1) * 8192) + j)
------array access p.addr (((i + 1) * 8192) + j)
------array access q.addr ((((i + 1) * 8192) + j) + 1)
----array access u.addr ((((i + 1) * 8192) + 8192) - 1)
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p.addr (((((i + 1) * 8192) + 8192) - 2) - j)
------array access u.addr (((((i + 1) * 8192) + 8192) - 1) - j)
------array access q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
------array access u.addr (((((i + 1) * 8192) + 8192) - 2) - j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 163
------Sample number: 26830
------Sample number: 26830
----Sample number: 163
------Sample number: 26830
------Sample number: 26830
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* v_addr0	98285 */
/* p_addr4	98285 */
/* u_addr19	98285 */
/* q_addr20	98285 */
/* u_addr6	98285 */
/* p_addr1	98285 */
/* v_addr2	98285 */
/* q_addr3	98285 */
/* q_addr11	98285 */
/* v_addr12	98285 */
/* p_addr13	98285 */
/* p_addr5	98285 */
/* v_addr16	98285 */
/* u_addr17	98285 */
/* p_addr18	98285 */
/* u_addr7	98285 */
/* u_addr8	98285 */
/* q_addr9	98285 */
/* p_addr10	98285 */
/* q_addr28	98285 */
/* u_addr29	98285 */
/* p_addr30	98285 */
/* u_addr31	98285 */
/* q_addr32	98285 */
/* v_addr14	98285 */
/* q_addr15	98285 */
/* p_addr21	98285 */
/* p_addr22	98285 */
/* v_addr23	98285 */
/* v_addr24	98285 */
/* v_addr25	98285 */
/* q_addr26	98285 */
/* p_addr27	98285 */
/* u_addr33	98285 */
#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <thread>
#include <unordered_map>
#include <vector>
#ifdef PAPI_TIMER
#  include "papi_timer.h"
#endif
using namespace std;
using namespace placeholders;
mutex mtx;
map<uint64_t, double> RT;
map<uint64_t, double> interceptRT;
map<uint64_t, double> scaleRT;
map<uint64_t, double> otherRT;
map<uint64_t, double> MR;
uint64_t sample_sum = 0;
double total_reuse = 0.0;
double total_src_neighbor = 0.0;
double total_sink_neighbor = 0.0;
double total_scale = 0.0;
double total_fold = 0.0;
double total_interchunk = 0.0;
double share_reuse = 0.0;
double total_smaller_reuse = 0.0;
struct Sample {
public:
    string name;
    vector<int> ivs;
    // int cid; // chunk this sample locates
    // int tid; // thread this sample belongs
    // int pos; // thread local position
    Sample() {}
    Sample(string ref, vector<int> iter) {
        name = ref;
        ivs = iter;
        // cid = floor(iter[0] / (CHUNK_SIZE * THREAD_NUM));
        // tid = iter[0] / CHUNK_SIZE - floor(iter[0] / (CHUNK_SIZE * THREAD_NUM)) * THREAD_NUM;
        // pos = iter[0] % CHUNK_SIZE;
    }
    string toString() {
        string s = "[" + name + "] {";
        vector<int>::iterator it;
        for (it = ivs.begin(); it != ivs.end(); ++it) {
            s += to_string(*it) + ",";
        }
        s.pop_back();
        s += "}";
        return s;
    }
    int compare(struct Sample other) {
        if (ivs.size() != other.ivs.size()) { return 2; }
        /* the same thread. compare the rest loop induction variables */ ;
        vector<int>::iterator selfit = ivs.begin();
        vector<int>::iterator otherit = other.ivs.begin();
        while (selfit != ivs.end() && otherit != other.ivs.end()) {
             if (*selfit < *otherit) {
                 return -1;
             } else if (*selfit > *otherit) {
                 return 1;
             }
             selfit++;
             otherit++;
        }
        /* all equal, these two samples are equal */ ;
        return 0;
    }
    bool operator==(const struct Sample & other) const {
    int i = 0;
    for (i = 0; i < ivs.size(); i++) {
        if (ivs[i] != other.ivs[i]) {
            return false;
        }
    }
        return name == other.name;
    }
};
// Hash function for Iteration
struct SampleHasher {
    size_t operator()(const struct Sample &iter) const {
        using std::string;
        using std::size_t;
        using std::hash;
       size_t hash_val = hash<string>()(iter.name);
        uint64_t bitmap = 0UL;
        int i = 2;
        for (auto iv : iter.ivs) {
            bitmap |= ((uint64_t)iv << (i * 14));
            i -= 1;
            if (i < 0) { break; }
        }
        hash_val ^= hash<uint64_t>()(bitmap);
        return hash_val;
    }
};
typedef struct Sample Sample;
void rtHistoCal( map<uint64_t, double> &rth, uint64_t rt, double val ) {
     unique_lock< mutex> lck (mtx, defer_lock);
    lck.lock();
    if (rth.find(rt) == rth.end()) { 
        rth[rt] = val;
    } else {
        rth[rt] += val;
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
bool isCompleteChunk(uint64_t is) {
    return ((is % (CHUNK_SIZE * THREAD_NUM)) == 0);
}
int getChunkNum(uint64_t is) {
    if (is % (CHUNK_SIZE * THREAD_NUM) != 0) {
        return is / (CHUNK_SIZE * THREAD_NUM) + 1;
    }
    return is / (CHUNK_SIZE * THREAD_NUM);
}
int getChunkID(uint64_t i) {
    return floor(i / (CHUNK_SIZE * THREAD_NUM));
}
int getThreadID(uint64_t i) {
    return i / CHUNK_SIZE - floor(i / (CHUNK_SIZE * THREAD_NUM))*THREAD_NUM ;
}
int getThreadLocalPos(uint64_t i) {
    return i % CHUNK_SIZE;
}
int search_src_candidate_neighbor(uint64_t i, function<uint64_t(uint64_t)> calAddr) {
    int c_start = i % CHUNK_SIZE + i / (CHUNK_SIZE * THREAD_NUM) * THREAD_NUM * CHUNK_SIZE;
    int c_end = c_start + (THREAD_NUM - 1) * CHUNK_SIZE;
    for (int c = i + CHUNK_SIZE; c <= c_end; c=c+CHUNK_SIZE) {
        if (calAddr(i) == calAddr(c)) { return getThreadID(c); }
    }
    return -1;
}
int search_sink_candidate_neighbor(uint64_t i, function<uint64_t(uint64_t)> calAddr) {
    int c_start = i % CHUNK_SIZE + i / (CHUNK_SIZE * THREAD_NUM) * THREAD_NUM * CHUNK_SIZE;
    for (int c = c_start; c <= i; c=c+CHUNK_SIZE) {
        if (calAddr(i) == calAddr(c)) { return getThreadID(c); }
    }
    return -1;
}
pair<uint64_t, int> parallel_predict(uint64_t i_src, uint64_t i_sink, uint64_t rt, uint64_t lsrc, uint64_t lsink, uint64_t middle_accesses, bool is_normal_ref, bool is_in_same_loop, int & type, function<uint64_t(uint64_t)> srcAddrCal, function<uint64_t(uint64_t)> sinkAddrCal) {
#ifdef DEBUG
    cout << "rt " << rt << endl;
#endif
    total_reuse += 1.0;
    uint64_t parallel_rt = rt;
    int tsrc = getThreadID(i_src);
    int tsink = getThreadID(i_sink);
    int dT = tsink - tsrc;
    int sink_neighbor_delta = 0;
    if (!is_in_same_loop || getChunkID(i_src) != getChunkID(i_sink)) {
#ifdef DEBUG
        cout << "Inter Chunk Reuse" << endl;
        ;
#endif
        type = 0; // code for inter chunk reuse
        total_interchunk += 1.0;
        if (dT != 0) { share_reuse += 1.0; }
        parallel_rt = rt * THREAD_NUM - CHUNK_SIZE * THREAD_NUM * (lsrc*(THREAD_NUM - tsrc) + lsink * tsink) + CHUNK_SIZE * THREAD_NUM * lsrc - (THREAD_NUM - 1) * middle_accesses + dT;
        if (parallel_rt < rt) { total_smaller_reuse += 1.0; }
        return make_pair(parallel_rt, 0);
    } else if (!is_normal_ref) {
        /* intra chunk reuse */
#ifdef DEBUG
        cout << "Neighboring Effect" << endl;
#endif
        int tsrc_neighbor = search_src_candidate_neighbor(i_src, srcAddrCal);
        int tsink_neighbor = search_sink_candidate_neighbor(i_sink, sinkAddrCal);
        if (tsrc_neighbor >= 0) {
#ifdef DEBUG
            cout << "Find sink in src neighbor at " << tsrc_neighbor << endl;
            cout << "Neighbor Effect: " << tsrc_neighbor - tsrc << endl;
#endif
            total_src_neighbor += 1.0;
            share_reuse += 1.0;
            type = 1; // code for src neighboring effect
        if ((tsrc_neighbor - tsrc) < rt) { total_smaller_reuse += 1.0; }
            return make_pair(tsrc_neighbor - tsrc, 1);
        } else if (tsink_neighbor >= 0) {
#ifdef DEBUG
            cout << "Find sink in sink neighbor at " << tsink_neighbor << endl;
            cout << "Neighbor Effect: " << rt * THREAD_NUM + tsink_neighbor - tsink << endl;
#endif
            if (getChunkID(i_src) == getChunkID(i_sink)) {
                total_sink_neighbor += 1.0;
                share_reuse += 1.0;
                type = 2; // code for sink neighboring effect
                if ((rt * THREAD_NUM + tsink_neighbor - tsink) < rt) { total_smaller_reuse += 1.0; }
                return make_pair(rt * THREAD_NUM + tsink_neighbor - tsink, 2);
            }
        }
    } else if (getChunkID(i_src) == getChunkID(i_sink)) {
        /* same thread -- scaling effect */
        if (dT == 0) {
#ifdef DEBUG
            cout << "Scaling Effect" << endl;
#endif
            if (sink_neighbor_delta == 0.0) { total_scale += 1.0; }
            parallel_rt = rt * THREAD_NUM; // * THREAD_NUM;
            type = 3; // code for scaling effect
        } else if (getThreadLocalPos(i_src) <= getThreadLocalPos(i_sink)) { // src-sink order
            if ((rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + dT) < 0) { printf("NORMAL ORDER NEGATIVE PRI\n"); }
#ifdef DEBUG
            cout << "Src-Sink Order Folding Effect" << endl;
#endif
            type = 4; // code for src-sink order folding effect
            if (sink_neighbor_delta == 0) {
                total_fold += 1.0;
                share_reuse += 1.0;
            }
            parallel_rt = rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + abs(dT);
            if (parallel_rt < rt) { total_smaller_reuse += 1.0; }
        } else { // sink-src order
            if ((rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + dT) < 0) { printf("REVERSE ORDER NEGATIVE PRI\n"); }
#ifdef DEBUG
            cout << "Sink-Src Order Folding Effect" << endl;
#endif
            // parallel_rt = CHUNK_SIZE * lsrc * THREAD_NUM * dT - (rt * THREAD_NUM) - abs(dT);
            type = 5; // code for sink-src order folding effect
            return make_pair(0UL, -1);
        }
    }
    return make_pair(parallel_rt + sink_neighbor_delta, type);
}
void RTtoMR_AET() {
     map<uint64_t, double> P;
    double total_num_RT = 0;
    uint64_t max_RT = 0;
    for ( map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
        total_num_RT += it->second;
        if (max_RT < it->first) {
            max_RT = it->first;
        }
    }
    double accumulate_num_RT = 0;
    for ( map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
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
void rtDistribution() {
    uint64_t maxInterceptRT = 0UL;
    int i = 0;
    map<uint64_t, double>::iterator mit = interceptRT.begin();
    for (; mit != interceptRT.end(); ++mit) {
        if (mit->first > maxInterceptRT) { maxInterceptRT = mit->first; }
    }
    for (; mit != interceptRT.end(); ++mit) {
        for (i = 1; i <= maxInterceptRT; i++) {
            rtHistoCal(RT, i, mit->second / maxInterceptRT);
        }
    }
    mit = scaleRT.begin();
    for (; mit != scaleRT.end(); ++mit) {
        uint64_t scaleRI = mit->first;
        for (i = (scaleRI / THREAD_NUM); i <= scaleRI; i++) {
            rtHistoCal(RT, i, mit->second / (scaleRI * (THREAD_NUM - 1) / THREAD_NUM));
        }
    }
    mit = otherRT.begin();
    for(; mit != otherRT.end(); ++mit) {
       rtHistoCal(RT, mit->first, mit->second);
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
     map<uint64_t, double>::iterator it1 = MR.begin();
     map<uint64_t, double>::iterator it2 = MR.begin();
    while(it1 != MR.end()) {
        while(1) {
             map<uint64_t, double>::iterator it3 = it2;
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
/* Dump the reuse statistics */
void statDump() {
    cout << "Total Neighboring (SRC) Reuses: " << total_src_neighbor / total_reuse << endl;
    cout << "Total Neighboring (SINK) Reuses: " << total_sink_neighbor / total_reuse << endl;
    cout << "Total Scaling Reuses: " << total_scale / total_reuse << endl;
    cout << "Total Folding Src-Sink Reuses: " << total_fold / total_reuse << endl;
    cout << "Total Inter Chunk Reuses: " << total_interchunk / total_reuse << endl;
    cout << "Total Intercept(Smaller) Reuses: " << total_smaller_reuse / total_reuse << endl;
    cout << "Total Shared Reuses: " << share_reuse / total_reuse << endl;
    cout << "Total Reuses: " << total_reuse << endl;
    return;
}
/* Array v_addr	i */ 
/* i */
/* v_addr ((0 + i) + 1) 0 */
int calAddrv_addr0( int i) {
    int result = (((0 + i) + 1)) * 8 / 64;
    return result;
}
/* Array p_addr	i */ 
/* i */
/* p_addr (((i + 1) * 8192) + 0) 1 */
int calAddrp_addr1( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array v_addr	i */ 
/* i */
/* v_addr (0 + (i + 1)) 2 */
int calAddrv_addr2( int i) {
    int result = ((0 + (i + 1))) * 8 / 64;
    return result;
}
/* Array q_addr	i */ 
/* i */
/* q_addr (((i + 1) * 8192) + 0) 3 */
int calAddrq_addr3( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + j) 4 */
int calAddrp_addr4( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr ((((i + 1) * 8192) + j) + 1) 5 */
int calAddrp_addr5( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array u_addr	j i */ 
/* i */
/* i */
/* i, j */
/* u_addr (((j + 1) * 8192) + i) 6 */
int calAddru_addr6( int i, int j) {
    int result = ((((j + 1) * 8192) + i)) * 8 / 64;
    return result;
}
/* Array u_addr	j i */ 
/* i */
/* i */
/* i, j */
/* u_addr ((((j + 1) * 8192) + i) + 1) 7 */
int calAddru_addr7( int i, int j) {
    int result = (((((j + 1) * 8192) + i) + 1)) * 8 / 64;
    return result;
}
/* Array u_addr	j i */ 
/* i */
/* i */
/* i, j */
/* u_addr ((((j + 1) * 8192) + i) + 2) 8 */
int calAddru_addr8( int i, int j) {
    int result = (((((j + 1) * 8192) + i) + 2)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr (((i + 1) * 8192) + j) 9 */
int calAddrq_addr9( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + j) 10 */
int calAddrp_addr10( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr ((((i + 1) * 8192) + j) + 1) 11 */
int calAddrq_addr11( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array v_addr	i */ 
/* i */
/* v_addr (67100672 + (i + 1)) 12 */
int calAddrv_addr12( int i) {
    int result = ((67100672 + (i + 1))) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + (8190 - j)) 13 */
int calAddrp_addr13( int i, int j) {
    int result = ((((i + 1) * 8192) + (8190 - j))) * 8 / 64;
    return result;
}
/* Array v_addr	j i */ 
/* i */
/* i */
/* i, j */
/* v_addr ((((8191 - j) * 8192) + i) + 1) 14 */
int calAddrv_addr14( int i, int j) {
    int result = (((((8191 - j) * 8192) + i) + 1)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr (((((i + 1) * 8192) + 8192) - 2) - j) 15 */
int calAddrq_addr15( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
/* Array v_addr	j i */ 
/* i */
/* i */
/* i, j */
/* v_addr (((8190 - j) * 8192) + (i + 1)) 16 */
int calAddrv_addr16( int i, int j) {
    int result = ((((8190 - j) * 8192) + (i + 1))) * 8 / 64;
    return result;
}
/* Array u_addr	i */ 
/* i */
/* u_addr (((i + 1) * 8192) + 0) 17 */
int calAddru_addr17( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array p_addr	i */ 
/* i */
/* p_addr (((i + 1) * 8192) + 0) 18 */
int calAddrp_addr18( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array u_addr	i */ 
/* i */
/* u_addr (((i + 1) * 8192) + 0) 19 */
int calAddru_addr19( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array q_addr	i */ 
/* i */
/* q_addr (((i + 1) * 8192) + 0) 20 */
int calAddrq_addr20( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + j) 21 */
int calAddrp_addr21( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr ((((i + 1) * 8192) + j) + 1) 22 */
int calAddrp_addr22( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array v_addr	i j */ 
/* i */
/* v_addr (((i * 8192) + j) + 1) 23 */
int calAddrv_addr23( int i, int j) {
    int result = ((((i * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array v_addr	i j */ 
/* i */
/* v_addr ((((i + 1) * 8192) + j) + 1) 24 */
int calAddrv_addr24( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array v_addr	i j */ 
/* i */
/* v_addr ((((i + 2) * 8192) + j) + 1) 25 */
int calAddrv_addr25( int i, int j) {
    int result = (((((i + 2) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr (((i + 1) * 8192) + j) 26 */
int calAddrq_addr26( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((i + 1) * 8192) + j) 27 */
int calAddrp_addr27( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr ((((i + 1) * 8192) + j) + 1) 28 */
int calAddrq_addr28( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* Array u_addr	i */ 
/* i */
/* u_addr ((((i + 1) * 8192) + 8192) - 1) 29 */
int calAddru_addr29( int i) {
    int result = (((((i + 1) * 8192) + 8192) - 1)) * 8 / 64;
    return result;
}
/* Array p_addr	i j */ 
/* i */
/* p_addr (((((i + 1) * 8192) + 8192) - 2) - j) 30 */
int calAddrp_addr30( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
/* Array u_addr	i j */ 
/* i */
/* u_addr (((((i + 1) * 8192) + 8192) - 1) - j) 31 */
int calAddru_addr31( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 1) - j)) * 8 / 64;
    return result;
}
/* Array q_addr	i j */ 
/* i */
/* q_addr (((((i + 1) * 8192) + 8192) - 2) - j) 32 */
int calAddrq_addr32( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
/* Array u_addr	i j */ 
/* i */
/* u_addr (((((i + 1) * 8192) + 8192) - 2) - j) 33 */
int calAddru_addr33( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
void ref_v_addr0() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("v_addr0", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr0( i);
                Sample iter("v_addr0", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr2( i);
                Sample iter("v_addr2", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
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
                int addr = calAddrv_addr12( i);
                Sample iter("v_addr12", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr14( i, j);
                    Sample iter("v_addr14", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr16( i, j);
                    Sample iter("v_addr16", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr23( i, j);
                    Sample iter("v_addr23", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr24( i, j);
                    Sample iter("v_addr24", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr25( i, j);
                    Sample iter("v_addr25", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr4() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr4", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr1( i);
                Sample iter("p_addr1", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr4( i, j);
                    Sample iter("p_addr4", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr5( i, j);
                    Sample iter("p_addr5", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr10( i, j);
                    Sample iter("p_addr10", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr13( i, j);
                    Sample iter("p_addr13", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
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
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_u_addr19() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("u_addr19", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr17( i);
                Sample iter("u_addr17", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr19( i);
                Sample iter("u_addr19", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
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
                int addr = calAddru_addr29( i);
                Sample iter("u_addr29", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr31( i, j);
                    Sample iter("u_addr31", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr33( i, j);
                    Sample iter("u_addr33", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr19, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
EndSample:
        return;
}
void ref_q_addr20() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr20", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr20( i);
                Sample iter("q_addr20", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr26( i, j);
                    Sample iter("q_addr26", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr28( i, j);
                    Sample iter("q_addr28", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr32( i, j);
                    Sample iter("q_addr32", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr20, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_u_addr6() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("u_addr6", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr6( i, j);
                    Sample iter("u_addr6", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr7( i, j);
                    Sample iter("u_addr7", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr8( i, j);
                    Sample iter("u_addr8", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr17( i);
                Sample iter("u_addr17", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr19( i);
                Sample iter("u_addr19", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
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
                int addr = calAddru_addr29( i);
                Sample iter("u_addr29", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr31( i, j);
                    Sample iter("u_addr31", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr33( i, j);
                    Sample iter("u_addr33", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr6, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr1() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr1", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr1( i);
                Sample iter("p_addr1", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr4( i, j);
                    Sample iter("p_addr4", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr5( i, j);
                    Sample iter("p_addr5", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr10( i, j);
                    Sample iter("p_addr10", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr13( i, j);
                    Sample iter("p_addr13", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
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
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_v_addr2() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("v_addr2", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr0( i);
                Sample iter("v_addr0", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr2( i);
                Sample iter("v_addr2", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
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
                int addr = calAddrv_addr12( i);
                Sample iter("v_addr12", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr14( i, j);
                    Sample iter("v_addr14", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr16( i, j);
                    Sample iter("v_addr16", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr23( i, j);
                    Sample iter("v_addr23", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr24( i, j);
                    Sample iter("v_addr24", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr25( i, j);
                    Sample iter("v_addr25", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr2, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_q_addr3() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr3", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr3( i);
                Sample iter("q_addr3", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr9( i, j);
                    Sample iter("q_addr9", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr11( i, j);
                    Sample iter("q_addr11", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr15( i, j);
                    Sample iter("q_addr15", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr20( i);
                Sample iter("q_addr20", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr26( i, j);
                    Sample iter("q_addr26", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr28( i, j);
                    Sample iter("q_addr28", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr32( i, j);
                    Sample iter("q_addr32", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr3, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_q_addr11() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr11", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr3( i);
                Sample iter("q_addr3", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr9( i, j);
                    Sample iter("q_addr9", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr11( i, j);
                    Sample iter("q_addr11", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr15( i, j);
                    Sample iter("q_addr15", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr20( i);
                Sample iter("q_addr20", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr26( i, j);
                    Sample iter("q_addr26", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr28( i, j);
                    Sample iter("q_addr28", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr32( i, j);
                    Sample iter("q_addr32", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr11, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_v_addr12() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("v_addr12", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr0( i);
                Sample iter("v_addr0", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr2( i);
                Sample iter("v_addr2", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
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
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr12( i);
                Sample iter("v_addr12", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr14( i, j);
                    Sample iter("v_addr14", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr16( i, j);
                    Sample iter("v_addr16", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr23( i, j);
                    Sample iter("v_addr23", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr24( i, j);
                    Sample iter("v_addr24", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr25( i, j);
                    Sample iter("v_addr25", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr12, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr13() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr13", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr1( i);
                Sample iter("p_addr1", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr4( i, j);
                    Sample iter("p_addr4", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr5( i, j);
                    Sample iter("p_addr5", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr10( i, j);
                    Sample iter("p_addr10", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr13( i, j);
                    Sample iter("p_addr13", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
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
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr13, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr5() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr5", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr1( i);
                Sample iter("p_addr1", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr4( i, j);
                    Sample iter("p_addr4", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr5( i, j);
                    Sample iter("p_addr5", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr10( i, j);
                    Sample iter("p_addr10", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr13( i, j);
                    Sample iter("p_addr13", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
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
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_v_addr16() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("v_addr16", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr0( i);
                Sample iter("v_addr0", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr2( i);
                Sample iter("v_addr2", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
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
                int addr = calAddrv_addr12( i);
                Sample iter("v_addr12", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr14( i, j);
                    Sample iter("v_addr14", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr16( i, j);
                    Sample iter("v_addr16", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr23( i, j);
                    Sample iter("v_addr23", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr24( i, j);
                    Sample iter("v_addr24", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr25( i, j);
                    Sample iter("v_addr25", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr16, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_u_addr17() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("u_addr17", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr17( i);
                Sample iter("u_addr17", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr19( i);
                Sample iter("u_addr19", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
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
                int addr = calAddru_addr29( i);
                Sample iter("u_addr29", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr31( i, j);
                    Sample iter("u_addr31", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr33( i, j);
                    Sample iter("u_addr33", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr17, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr18() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr18", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr18, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_u_addr7() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("u_addr7", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr6( i, j);
                    Sample iter("u_addr6", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr7( i, j);
                    Sample iter("u_addr7", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr8( i, j);
                    Sample iter("u_addr8", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr17( i);
                Sample iter("u_addr17", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr19( i);
                Sample iter("u_addr19", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
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
                int addr = calAddru_addr29( i);
                Sample iter("u_addr29", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr31( i, j);
                    Sample iter("u_addr31", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr33( i, j);
                    Sample iter("u_addr33", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr7, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
EndSample:
        return;
}
void ref_u_addr8() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("u_addr8", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr6( i, j);
                    Sample iter("u_addr6", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr6, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr7( i, j);
                    Sample iter("u_addr7", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr7, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr8( i, j);
                    Sample iter("u_addr8", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr8, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr17( i);
                Sample iter("u_addr17", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr19( i);
                Sample iter("u_addr19", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
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
                int addr = calAddru_addr29( i);
                Sample iter("u_addr29", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr31( i, j);
                    Sample iter("u_addr31", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr33( i, j);
                    Sample iter("u_addr33", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr8, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
EndSample:
        return;
}
void ref_q_addr9() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr9", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr3( i);
                Sample iter("q_addr3", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr9( i, j);
                    Sample iter("q_addr9", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr11( i, j);
                    Sample iter("q_addr11", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr15( i, j);
                    Sample iter("q_addr15", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr20( i);
                Sample iter("q_addr20", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr26( i, j);
                    Sample iter("q_addr26", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr28( i, j);
                    Sample iter("q_addr28", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr32( i, j);
                    Sample iter("q_addr32", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr9, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr10() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr10", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr1( i);
                Sample iter("p_addr1", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr1, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr4( i, j);
                    Sample iter("p_addr4", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr5( i, j);
                    Sample iter("p_addr5", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr5, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr10( i, j);
                    Sample iter("p_addr10", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr10, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr13( i, j);
                    Sample iter("p_addr13", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr13, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
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
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr10, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_q_addr28() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr28", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr20( i);
                Sample iter("q_addr20", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr26( i, j);
                    Sample iter("q_addr26", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr28( i, j);
                    Sample iter("q_addr28", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr32( i, j);
                    Sample iter("q_addr32", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr28, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_u_addr29() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 163;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("u_addr29", {iSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr17( i);
                Sample iter("u_addr17", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr19( i);
                Sample iter("u_addr19", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
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
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr29( i);
                Sample iter("u_addr29", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr31( i, j);
                    Sample iter("u_addr31", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr33( i, j);
                    Sample iter("u_addr33", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr29, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr30() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr30", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr30, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_u_addr31() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("u_addr31", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr17( i);
                Sample iter("u_addr17", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr19( i);
                Sample iter("u_addr19", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
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
                int addr = calAddru_addr29( i);
                Sample iter("u_addr29", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr31( i, j);
                    Sample iter("u_addr31", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr33( i, j);
                    Sample iter("u_addr33", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr31, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
EndSample:
        return;
}
void ref_q_addr32() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr32", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr20( i);
                Sample iter("q_addr20", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr26( i, j);
                    Sample iter("q_addr26", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr28( i, j);
                    Sample iter("q_addr28", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr32( i, j);
                    Sample iter("q_addr32", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr32, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_v_addr14() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("v_addr14", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr0( i);
                Sample iter("v_addr0", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr0, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrv_addr2( i);
                Sample iter("v_addr2", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr2, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
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
                int addr = calAddrv_addr12( i);
                Sample iter("v_addr12", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr12, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr14( i, j);
                    Sample iter("v_addr14", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr14, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr16( i, j);
                    Sample iter("v_addr16", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr16, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr23( i, j);
                    Sample iter("v_addr23", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr24( i, j);
                    Sample iter("v_addr24", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr25( i, j);
                    Sample iter("v_addr25", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: false */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr14, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_q_addr15() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr15", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr3( i);
                Sample iter("q_addr3", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr3, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr9( i, j);
                    Sample iter("q_addr9", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr9, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr11( i, j);
                    Sample iter("q_addr11", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr11, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr15( i, j);
                    Sample iter("q_addr15", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr15, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr20( i);
                Sample iter("q_addr20", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: false */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        if (isCompleteChunk(8190)) {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        } else {
                                    middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                        }
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr26( i, j);
                    Sample iter("q_addr26", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr28( i, j);
                    Sample iter("q_addr28", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr32( i, j);
                    Sample iter("q_addr32", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: false */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr15, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8190)) {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            } else {
                                            middle_accesses += ((getChunkNum(8190) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8190 % (THREAD_NUM * CHUNK_SIZE))) * 98285 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 98285;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, false, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr21() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr21", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr21, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr22() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr22", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr22, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_v_addr23() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("v_addr23", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr23( i, j);
                    Sample iter("v_addr23", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr24( i, j);
                    Sample iter("v_addr24", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr25( i, j);
                    Sample iter("v_addr25", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr23, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_v_addr24() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("v_addr24", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr23( i, j);
                    Sample iter("v_addr23", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr24( i, j);
                    Sample iter("v_addr24", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr25( i, j);
                    Sample iter("v_addr25", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr24, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_v_addr25() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("v_addr25", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr23( i, j);
                    Sample iter("v_addr23", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr23, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr24( i, j);
                    Sample iter("v_addr24", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr24, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrv_addr25( i, j);
                    Sample iter("v_addr25", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrv_addr25, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrv_addr25, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_q_addr26() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr26", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr20( i);
                Sample iter("q_addr20", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr20, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr26( i, j);
                    Sample iter("q_addr26", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr26, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr28( i, j);
                    Sample iter("q_addr28", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr28, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr32( i, j);
                    Sample iter("q_addr32", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr26, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr32, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_p_addr27() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("p_addr27", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrp_addr18( i);
                Sample iter("p_addr18", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr18, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr21( i, j);
                    Sample iter("p_addr21", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr21, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr22( i, j);
                    Sample iter("p_addr22", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr22, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr27( i, j);
                    Sample iter("p_addr27", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr27, _1, j);
                            /* i 98285 */
                            /* j 8 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
            }
            }
            if (cntStart == true) cnt++;
            {
            int jLB5 = 0;
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrp_addr30( i, j);
                    Sample iter("p_addr30", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr27, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr30, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
            }
        }
        }
EndSample:
        return;
}
void ref_u_addr33() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26830;) {
SAMPLE:
        int iSample = rand() % (8190 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8190) { goto SAMPLE; }
        if ( (8190 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8190 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("u_addr33", {iSample, jSample});
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
        cout << "Samples:" << sample.toString() << endl;
#endif
        samples[sample] = 1;
        if (s == 0 || sample.compare(sStart) < 0) { sStart = sample; }
        s += 1UL;
        }
        int i_Start = sStart.ivs[0];
        int j_Start = sStart.ivs[1];
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB3 = i_Start;
        for ( int i = iLB3; i < 8190; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr17( i);
                Sample iter("u_addr17", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr17, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                int addr = calAddru_addr19( i);
                Sample iter("u_addr19", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr19, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            if (cntStart == true) cnt++;
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8190; j=(j + 1)) {
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
                int addr = calAddru_addr29( i);
                Sample iter("u_addr29", {i});
                if (LAT.find(addr) != LAT.end()) {
                        /* Find a reuse. Find the Sample object of the reuse src */
                        uint64_t ris = cnt - LAT[addr];
                        Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                        cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                        /* is_src_loop_outermost: 1 */
                        /* is_sink_loop_outermost: 1 */
                        /* is_normal_ref: false */
                        /* is_in_same_loop: true */
                        /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr29, _1);
                        /* i 98285 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
                        rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                    LAT.erase(addr);
                    samples.erase(srcSample);
                    if (samples.size() == 0) { goto EndSample; }
                } // end of check addr in LAT
                if (samples.find(iter) != samples.end()) {
                    LAT[addr] = cnt;
                    LATSampleIterMap[cnt] = iter;
                }
            }
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 8190; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr31( i, j);
                    Sample iter("u_addr31", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr31, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddru_addr33( i, j);
                    Sample iter("u_addr33", {i, j});
                    if (LAT.find(addr) != LAT.end()) {
                            /* Find a reuse. Find the Sample object of the reuse src */
                            uint64_t ris = cnt - LAT[addr];
                            Sample srcSample = LATSampleIterMap[LAT[addr]];
#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)
                            cout << "[" << cnt - LAT[addr] << "] " <<  srcSample.toString() << " -> " << iter.toString() << endl;
#endif
                            /* is_src_loop_outermost: 1 */
                            /* is_sink_loop_outermost: 1 */
                            /* is_normal_ref: true */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddru_addr33, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddru_addr33, _1, j);
                            /* i 98285 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 98285, 98285, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
                            rtHistoCal(RT, get<0>(parallel_rt), 1.0);
                        LAT.erase(addr);
                        samples.erase(srcSample);
                        if (samples.size() == 0) { goto EndSample; }
                    } // end of check addr in LAT
                    if (samples.find(iter) != samples.end()) {
                        LAT[addr] = cnt;
                        LATSampleIterMap[cnt] = iter;
                    }
                }
            }
            }
        }
        }
EndSample:
        return;
}
int main() {
#ifdef PAPI_TIMER
    PAPI_timer_init();
    PAPI_timer_start();
#endif
    thread t_v_addr_0(ref_v_addr0);
    thread t_p_addr_4(ref_p_addr4);
    thread t_u_addr_19(ref_u_addr19);
    thread t_q_addr_20(ref_q_addr20);
    thread t_u_addr_6(ref_u_addr6);
    thread t_p_addr_1(ref_p_addr1);
    thread t_v_addr_2(ref_v_addr2);
    thread t_q_addr_3(ref_q_addr3);
    thread t_q_addr_11(ref_q_addr11);
    thread t_v_addr_12(ref_v_addr12);
    thread t_p_addr_13(ref_p_addr13);
    thread t_p_addr_5(ref_p_addr5);
    thread t_v_addr_16(ref_v_addr16);
    thread t_u_addr_17(ref_u_addr17);
    thread t_p_addr_18(ref_p_addr18);
    thread t_u_addr_7(ref_u_addr7);
    thread t_u_addr_8(ref_u_addr8);
    thread t_q_addr_9(ref_q_addr9);
    thread t_p_addr_10(ref_p_addr10);
    thread t_q_addr_28(ref_q_addr28);
    thread t_u_addr_29(ref_u_addr29);
    thread t_p_addr_30(ref_p_addr30);
    thread t_u_addr_31(ref_u_addr31);
    thread t_q_addr_32(ref_q_addr32);
    thread t_v_addr_14(ref_v_addr14);
    thread t_q_addr_15(ref_q_addr15);
    thread t_p_addr_21(ref_p_addr21);
    thread t_p_addr_22(ref_p_addr22);
    thread t_v_addr_23(ref_v_addr23);
    thread t_v_addr_24(ref_v_addr24);
    thread t_v_addr_25(ref_v_addr25);
    thread t_q_addr_26(ref_q_addr26);
    thread t_p_addr_27(ref_p_addr27);
    thread t_u_addr_33(ref_u_addr33);
    t_v_addr_0.join();
    t_p_addr_4.join();
    t_u_addr_19.join();
    t_q_addr_20.join();
    t_u_addr_6.join();
    t_p_addr_1.join();
    t_v_addr_2.join();
    t_q_addr_3.join();
    t_q_addr_11.join();
    t_v_addr_12.join();
    t_p_addr_13.join();
    t_p_addr_5.join();
    t_v_addr_16.join();
    t_u_addr_17.join();
    t_p_addr_18.join();
    t_u_addr_7.join();
    t_u_addr_8.join();
    t_q_addr_9.join();
    t_p_addr_10.join();
    t_q_addr_28.join();
    t_u_addr_29.join();
    t_p_addr_30.join();
    t_u_addr_31.join();
    t_q_addr_32.join();
    t_v_addr_14.join();
    t_q_addr_15.join();
    t_p_addr_21.join();
    t_p_addr_22.join();
    t_v_addr_23.join();
    t_v_addr_24.join();
    t_v_addr_25.join();
    t_q_addr_26.join();
    t_p_addr_27.join();
    t_u_addr_33.join();
    RTtoMR_AET();
#ifdef PAPI_TIMER
    PAPI_timer_end();
    PAPI_timer_print();
#endif
    statDump();
    cout << "Samples: " << "1609908300" << endl;
    rtDump();
    dumpMR();
    return 0;
}
 /* Analyze function: adi */ 
