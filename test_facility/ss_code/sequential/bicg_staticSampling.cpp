
 /* Start to analysis array index
Array index info: Total number of references: 10
s.addr i
q.addr i
q.addr i
A.addr ((i * 8192) + j)
p.addr j
q.addr i
s.addr i
r.addr j
A.addr ((j * 8192) + i)
s.addr i
BC Array cost info: Total number of arrays: 5
A.addr 4
p.addr 2
q.addr 8
r.addr 2
s.addr 8
BC Average cost: 4.800000e+00

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
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----array access s.addr i
--i
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----array access q.addr i
----j
----Loop Bound: (0, 8192)
----Loop inc: (j + 1)
----Loop predicate: <
------array access q.addr i
------array access A.addr ((i * 8192) + j)
------array access p.addr j
------array access q.addr i
--i
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 8192)
----Loop inc: (j + 1)
----Loop predicate: <
------array access s.addr i
------array access r.addr j
------array access A.addr ((j * 8192) + i)
------array access s.addr i

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 163
----Sample number: 163
------Sample number: 26843
----Sample number: 163
------Sample number: 26843
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* s_addr0	1 */
/* q_addr1	32769 */
/* q_addr2	32769 */
/* A_addr3	32769 */
/* p_addr4	32769 */
/* q_addr5	32769 */
/* s_addr6	32768 */
/* r_addr7	32768 */
/* A_addr8	32768 */
/* s_addr9	32768 */
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
/* Array s_addr	i */ 
/* i */
/* s_addr i 0 */
int calAddrs_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array q_addr	i */ 
/* i */
/* q_addr i 1 */
int calAddrq_addr1( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array q_addr	i */ 
/* i */
/* q_addr i 2 */
int calAddrq_addr2( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array A_addr	i j */ 
/* i */
/* A_addr ((i * 8192) + j) 3 */
int calAddrA_addr3( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* Array p_addr	j */ 
/* i */
/* p_addr j 4 */
int calAddrp_addr4( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* Array q_addr	i */ 
/* i */
/* q_addr i 5 */
int calAddrq_addr5( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array s_addr	i */ 
/* i */
/* s_addr i 6 */
int calAddrs_addr6( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* Array r_addr	j */ 
/* i */
/* r_addr j 7 */
int calAddrr_addr7( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* Array A_addr	j i */ 
/* i */
/* i */
/* i, j */
/* A_addr ((j * 8192) + i) 8 */
int calAddrA_addr8( int i, int j) {
    int result = (((j * 8192) + i)) * 8 / 64;
    return result;
}
/* Array s_addr	i */ 
/* i */
/* s_addr i 9 */
int calAddrs_addr9( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
void ref_s_addr0() {
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
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("s_addr0", {iSample});
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
        for ( int i = iLB0; i < 8192; i=(i + 1)) {
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrs_addr0( i);
                Sample iter("s_addr0", {i});
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
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrs_addr0, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrs_addr0, _1);
                        /* i 1 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 1, 1, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
        {
        int iLB1 = 0;
        for ( int i = iLB1; i < 8192; i=(i + 1)) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            for ( int j = jLB2; j < 8192; j=(j + 1)) {
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
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrs_addr6( i, j);
                    Sample iter("s_addr6", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrs_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrs_addr6, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192)) {
                                            middle_accesses += ((getChunkNum(8192) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 1 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8192 % (THREAD_NUM * CHUNK_SIZE))) * 1 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 32768;
                            }
                            middle_accesses += 268443648;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 1, 32768, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
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
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrs_addr9( i, j);
                    Sample iter("s_addr9", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrs_addr0, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrs_addr9, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192)) {
                                            middle_accesses += ((getChunkNum(8192) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 1 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8192 % (THREAD_NUM * CHUNK_SIZE))) * 1 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 32768;
                            }
                            middle_accesses += 268443648;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 1, 32768, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
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
void ref_q_addr1() {
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
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        string idx_string =  to_string(iSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr1", {iSample});
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 8192; i=(i + 1)) {
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr1( i);
                Sample iter("q_addr1", {i});
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
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr1, _1);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr1, _1);
                        /* i 32769 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
            for ( int j = jLB2; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr2( i, j);
                    Sample iter("q_addr2", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr2, _1, j);
                            /* i 32769 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr5( i, j);
                    Sample iter("q_addr5", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr1, _1);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr5, _1, j);
                            /* i 32769 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
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
void ref_q_addr2() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26843;) {
SAMPLE:
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8192 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr2", {iSample, jSample});
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 8192; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr1( i);
                Sample iter("q_addr1", {i});
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
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr2, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr1, _1);
                        /* i 32769 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
            for ( int j = jLB2; j < 8192; j=(j + 1)) {
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr2( i, j);
                    Sample iter("q_addr2", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr2, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr2, _1, j);
                            /* i 32769 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr5( i, j);
                    Sample iter("q_addr5", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr2, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr5, _1, j);
                            /* i 32769 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
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
void ref_A_addr3() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26843;) {
SAMPLE:
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8192 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("A_addr3", {iSample, jSample});
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 8192; i=(i + 1)) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrA_addr3( i, j);
                    Sample iter("A_addr3", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr3, _1, j);
                            /* i 32769 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, true, true, reuse_type, srcAddrCal, sinkAddrCal);
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
            }
            }
        }
        }
        {
        int iLB3 = 0;
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrA_addr8( i, j);
                    Sample iter("A_addr8", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr3, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            if (isCompleteChunk(8192)) {
                                            middle_accesses += ((getChunkNum(8192) - getChunkID((srcSample.ivs[0] - 0)) - 1) * (CHUNK_SIZE * THREAD_NUM)) * 32769 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 32768;
                            } else {
                                            middle_accesses += ((getChunkNum(8192) - getChunkID((srcSample.ivs[0] - 0)) - 2) * CHUNK_SIZE * THREAD_NUM + (8192 % (THREAD_NUM * CHUNK_SIZE))) * 32769 + getChunkID((i - 0)) * CHUNK_SIZE * THREAD_NUM * 32768;
                            }
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32768, middle_accesses, false, false, reuse_type, srcAddrCal, sinkAddrCal);
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
void ref_p_addr4() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26843;) {
SAMPLE:
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8192 - 0) + 0;
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 8192; i=(i + 1)) {
            if (cntStart == true) cnt++;
            {
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
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
                            /* is_normal_ref: false */
                            /* is_in_same_loop: true */
                            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrp_addr4, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrp_addr4, _1, j);
                            /* i 32769 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
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
void ref_q_addr5() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26843;) {
SAMPLE:
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8192 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("q_addr5", {iSample, jSample});
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 8192; i=(i + 1)) {
            if (cntStart == true) {
                cnt++;
                int addr = calAddrq_addr1( i);
                Sample iter("q_addr1", {i});
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
                        function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr5, _1, srcSample.ivs[1]);
                        function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr1, _1);
                        /* i 32769 */
                        /* compute the number of accesses between source and sink chunk */
                        uint64_t middle_accesses = 0;
                        middle_accesses += 0;
#ifdef DEBUG
                        cout << " middle_access is " << middle_accesses << endl;
#endif
                        int reuse_type = -1;
                        pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
            for ( int j = jLB2; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr2( i, j);
                    Sample iter("q_addr2", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr2, _1, j);
                            /* i 32769 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrq_addr5( i, j);
                    Sample iter("q_addr5", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrq_addr5, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrq_addr5, _1, j);
                            /* i 32769 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32769, 32769, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
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
void ref_s_addr6() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26843;) {
SAMPLE:
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8192 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("s_addr6", {iSample, jSample});
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
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrs_addr6( i, j);
                    Sample iter("s_addr6", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrs_addr6, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrs_addr6, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32768, 32768, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrs_addr9( i, j);
                    Sample iter("s_addr9", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrs_addr6, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrs_addr9, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32768, 32768, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
void ref_r_addr7() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26843;) {
SAMPLE:
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8192 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("r_addr7", {iSample, jSample});
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
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrr_addr7( i, j);
                    Sample iter("r_addr7", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrr_addr7, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrr_addr7, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32768, 32768, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
            }
            }
        }
        }
EndSample:
        return;
}
void ref_A_addr8() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26843;) {
SAMPLE:
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8192 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("A_addr8", {iSample, jSample});
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
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrA_addr8( i, j);
                    Sample iter("A_addr8", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrA_addr8, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrA_addr8, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32768, 32768, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
void ref_s_addr9() {
    /* Generating sampling loop */
    set<string> record;
    // access time -> Sample Object
    unordered_map<uint64_t, Sample> LATSampleIterMap;
    // address -> access time
    unordered_map<int, uint64_t> LAT;
    unordered_map<Sample, int, SampleHasher> samples;
    Sample sStart;
    for ( int s = 0; s < 26843;) {
SAMPLE:
        int iSample = rand() % (8192 - 0) + 0;
        if (iSample % 1 != 0) goto SAMPLE; 
        if (iSample + THREAD_NUM * CHUNK_SIZE > 8192) { goto SAMPLE; }
        if ( (8192 - 0) == 0) goto SAMPLE;
        int jSample = rand() % (8192 - 0) + 0;
        if (jSample % 1 != 0) goto SAMPLE; 
        string idx_string =  to_string(iSample) + "_" +  to_string(jSample) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        Sample sample("s_addr9", {iSample, jSample});
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
        for ( int i = iLB3; i < 8192; i=(i + 1)) {
            {
            int jLB4 = 0;
            if ( i == i_Start ) {
                jLB4 = j_Start;
            }
            for ( int j = jLB4; j < 8192; j=(j + 1)) {
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrs_addr6( i, j);
                    Sample iter("s_addr6", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrs_addr9, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrs_addr6, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32768, 32768, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    int addr = calAddrs_addr9( i, j);
                    Sample iter("s_addr9", {i, j});
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
                            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddrs_addr9, _1, srcSample.ivs[1]);
                            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddrs_addr9, _1, j);
                            /* i 32768 */
                            /* j 4 */
                            /* compute the number of accesses between source and sink chunk */
                            uint64_t middle_accesses = 0;
                            middle_accesses += 0;
#ifdef DEBUG
                            cout << " middle_access is " << middle_accesses << endl;
#endif
                            int reuse_type = -1;
                            pair<uint64_t, int> parallel_rt = parallel_predict((srcSample.ivs[0] - 0), (i - 0), ris, 32768, 32768, middle_accesses, false, true, reuse_type, srcAddrCal, sinkAddrCal);
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
    thread t_s_addr_0(ref_s_addr0);
    thread t_q_addr_1(ref_q_addr1);
    thread t_q_addr_2(ref_q_addr2);
    thread t_A_addr_3(ref_A_addr3);
    thread t_p_addr_4(ref_p_addr4);
    thread t_q_addr_5(ref_q_addr5);
    thread t_s_addr_6(ref_s_addr6);
    thread t_r_addr_7(ref_r_addr7);
    thread t_A_addr_8(ref_A_addr8);
    thread t_s_addr_9(ref_s_addr9);
    t_s_addr_0.join();
    t_q_addr_1.join();
    t_q_addr_2.join();
    t_A_addr_3.join();
    t_p_addr_4.join();
    t_q_addr_5.join();
    t_s_addr_6.join();
    t_r_addr_7.join();
    t_A_addr_8.join();
    t_s_addr_9.join();
    RTtoMR_AET();
#ifdef PAPI_TIMER
    PAPI_timer_end();
    PAPI_timer_print();
#endif
    statDump();
    cout << "Samples: " << "536887296" << endl;
    rtDump();
    dumpMR();
    return 0;
}
 /* Analyze function: bicg_cpu */ 
