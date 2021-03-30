#include "uiAccCodeGen_ref.hpp"


namespace uiAccCodeGen_ref {
    char AccLevelUISamplingCodeGen_ref::ID = 0;
    static RegisterPass<AccLevelUISamplingCodeGen_ref> X("uiAccCodeGen_ref", "parallel 1:1 uniform interleaving static sampling code generating pass (access based)", false, false);

    AccLevelUISamplingCodeGen_ref::AccLevelUISamplingCodeGen_ref() : FunctionPass(ID) {}

    void AccLevelUISamplingCodeGen_ref::numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree->AA != NULL) {
            if (refToSameArrayCnt.find(arrayName[LoopRefTree->AA]) == refToSameArrayCnt.end()) {
                refNumber[LoopRefTree->AA] = 0;
                refToSameArrayCnt[arrayName[LoopRefTree->AA]] = 1;
            } else {
                refNumber[LoopRefTree->AA] = refToSameArrayCnt[arrayName[LoopRefTree->AA]];
                refToSameArrayCnt[arrayName[LoopRefTree->AA]] += 1;
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                numberRefToSameArray(*it);
            }
        }
        
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree->L != NULL) {
            if (loopNumber.find(LoopRefTree->L) == loopNumber.end()) {
                loopNumber[LoopRefTree->L] = loopNumber.size();
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                if ((*it)->isThreadNode) { continue; }
                numberLoops(*it);
            }
        }
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree == NULL) {
            return;
        }
        
        if (LoopRefTree->L != NULL) {
            for (std::vector<Value *>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvName[*it] = (*it)->getName().str();
                if (std::find(indvName[*it].begin(), indvName[*it].end(), '.') != indvName[*it].end()) {
                    indvName[*it].erase(std::find(indvName[*it].begin(), indvName[*it].end(), '.'));
                }
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                initIndvName(*it);
            }
        }
        
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::initArrayName() {
        
        for (std::map<Instruction*, std::string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            it->second.replace(std::find(it->second.begin(), it->second.end(), '.'), std::find(it->second.begin(), it->second.end(), '.') +1, 1, '_');
        }
        
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<string> indvs) {
        
        if (LoopRefTree->L != NULL) {
            for (std::vector<Value*>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvs.push_back(indvName[(*it)]);
            }
        }
        
        if (LoopRefTree->AA != NULL) {
			errs() << "/* " + arrayName[LoopRefTree->AA] + " " + arrayExpression[LoopRefTree->AA] + " " + std::to_string(refNumber[LoopRefTree->AA])  + " */\n";
            errs() << "int calAddr" + arrayName[LoopRefTree->AA] + std::to_string(refNumber[LoopRefTree->AA]) + "( ";
            string arguments = "";
            for (std::vector<string>::iterator it = indvs.begin(), eit = indvs.end(); it != eit; ++it) {
                arguments += "int " + (*it) + ", ";
            }
            if (arguments.size() != 0) {
                arguments.pop_back();
                arguments.pop_back();
            }
            errs() << arguments;
            errs() << ") {\n";
#if defined (CLS) && defined (DS)
            errs() << "    int result = (" + arrayExpression[LoopRefTree->AA] + ") * " + std::to_string(DS) + " / " + std::to_string(CLS) + ";\n";
#else
            errs() << "    int result = " + arrayExpression[LoopRefTree->AA] + ";\n";
#endif
            errs() << "    return result;\n";
            errs() << "}\n";
            return;
        }

        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                addrCalFuncGen(*it, indvs);
            }
        }
        
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        std::vector<string> indvs;
        addrCalFuncGen(LoopRefTree, indvs);
        
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::headerGen() {
        
        errs() << "#include <map>\n";
        errs() << "#include <set>\n";
        errs() << "#include <vector>\n";
        errs() << "#include <cstdlib>\n";
        errs() << "#include <iostream>\n";
        errs() << "#include <cmath>\n";
        errs() << "#include <numeric>\n";
        errs() << "#include <utility>\n";
#ifdef PARALLEL_CXX_THREAD
        errs() << "#include <thread>\n";
        errs() << "#include <mutex>\n";
#endif
#ifdef REFERENCE_GROUP
        errs() << "#include <algorithm>\n";
#endif
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "#include \"papi_timer.h\"\n";
        errs() << "#endif\n";
#ifdef PARALLEL
        errs() << "#ifndef THREAD_NUM\n";
        errs() << "#    define THREAD_NUM   4\n";
        errs() << "#endif\n";
        errs() << "#ifndef BIN_SIZE\n";
        errs() << "#    define BIN_SIZE   4\n";
        errs() << "#endif\n";
        errs() << "#ifndef CHUNK_SIZE\n";
        errs() << "#    define CHUNK_SIZE   4\n";
        errs() << "#endif\n";
#endif

        errs() << "using namespace std;\n";
        
#ifdef PARALLEL_CXX_THREAD
        errs() << "std::mutex mtx;\n";
#endif
#ifdef REFERENCE_GROUP
        errs() << "std::map<string, map<uint64_t, double>> refRT;\n";
#endif
        errs() << "std::map<uint64_t, double> RT;\n";
        errs() << "std::map<uint64_t, double> MR;\n";
        return;
    }

#ifdef DumpRTMR
    void AccLevelUISamplingCodeGen_ref::rtHistoGen() {
        errs() << "void rtHistoCal( map<uint64_t, double> &rth, uint64_t rt, double val ) {\n";
        errs() << "    if ( val <= 0) {\n;";
        errs() << "        return;\n";
        errs() << "    }\n";
#ifdef PARALLEL_OMP
        errs() << "    #pragma omp critical\n";
        errs() << "    {";
        errs() << "        if (rth.find(rt) == rth.end()) { \n";
        errs() << "            rth[rt] = val;\n";
        errs() << "        } else {\n";
        errs() << "            rth[rt] += val;\n";
        errs() << "        }\n";
        errs() << "    }";
#elif defined(PARALLEL_CXX_THREAD)
        errs() << "    std::unique_lock<std::mutex> lck (mtx,std::defer_lock);\n";
        errs() << "    lck.lock();\n";
#endif
        errs() << "    if (rth.find(rt) == rth.end()) { \n";
        errs() << "        rth[rt] = val;\n";
        errs() << "    } else {\n";
        errs() << "        rth[rt] += val;\n";
        errs() << "    }\n";
#if defined(PARALLEL_CXX_THREAD)
        errs() << "    lck.unlock();\n";
#endif
        errs() << "    return;\n";
        errs() << "}\n";
#ifdef REFERENCE_GROUP   
        errs() << "\n";    
        errs() << "void refRTHistoCal(map<string, map<uint64_t, double>> &rth, uint64_t rt, double val, string ref ) {\n";
        errs() << "    if ( val <= 0) {\n;";
        errs() << "        return;\n";
        errs() << "    }\n";
#ifdef PARALLEL_CXX_THREAD
        errs() << "    unique_lock< mutex> lck (mtx, defer_lock);\n";
        errs() << "    lck.lock();\n";
#endif
        errs() << "    if (rth.find(ref) == rth.end()) { \n";
        errs() << "        rth[ref][rt] = val;\n";
        errs() << "    }\n";
        errs() << "    else {\n";
        errs() << "        if (rth[ref].find(rt) == rth[ref].end()) { \n";
        errs() << "            rth[ref][rt] = val;\n";
        errs() << "        } else {\n";
        errs() << "            rth[ref][rt] += val;\n";
        errs() << "        }\n";
        errs() << "    }\n";
#ifdef PARALLEL_CXX_THREAD
        errs() << "    lck.unlock();\n";
#endif
        errs() << "    return;\n";
        errs() << "}\n";
#endif
    }

#elif defined(DumpRefLease)
    void AccLevelUISamplingCodeGen_ref::rtHistoGen() {
        errs() << "map<uint64_t, map<uint64_t, uint64_t>* > RI;\n";
        errs() << "map<uint64_t, map<uint64_t, double>* > hits;\n";
        errs() << "map<uint64_t, map<uint64_t, double>* > costs;\n";
        errs() << "map<uint64_t, double> sampledCnt;\n";
        errs() << "map<uint64_t, double> accessRatio;\n";
        errs() << "map<uint64_t, uint64_t> Lease;\n";

        errs() << "void rtHistoCal(uint64_t ri, uint64_t ref_id) {\n";
        errs() << "    if (RI.find(ref_id) != RI.end()) {\n";
        errs() << "        if ((*RI[ref_id]).find(ri) != (*RI[ref_id]).end()) {\n";
        errs() << "            (*RI[ref_id])[ri] ++;\n";
        errs() << "        } else {\n";
        errs() << "            (*RI[ref_id])[ri] = 1;\n";
        errs() << "        }\n";
        errs() << "    } else {\n";
        errs() << "        RI[ref_id] = new map<uint64_t, uint64_t>;\n";
        errs() << "        (*RI[ref_id])[ri] = 1;\n";
        errs() << "    }\n";
        errs() << "\n";
        errs() << "    // Init leases to all references to be 0\n";
        errs() << "    if (Lease.find(ref_id) == Lease.end()) {\n";
        errs() << "        Lease[ref_id] = 0;\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
    }
#endif

    /* Generate the function to calculate the bins */
    void AccLevelUISamplingCodeGen_ref::subBlkRTGen() {
        string space = "    ";
        errs() << "void subBlkRT(map<uint64_t, double> &rth, int rt, double cnt) {\n";
        errs() << space + "int msb = 0;\n";
        errs() << space + "int tmp_rt = rt;\n";
        errs() << space + "while(tmp_rt != 0) {\n";
        errs() << space + "    tmp_rt = tmp_rt / 2;\n";
        errs() << space + "    ++msb;\n";
        errs() << space + "}\n";
        errs() << space + "if (msb >= BIN_SIZE) {\n";
        errs() << space + "    int diff = (pow(2, msb) - pow(2, msb-1)) / BIN_SIZE;\n";
        errs() << space + "    for (int b = pow(2, msb-1); b <= pow(2, msb); b+=diff) {\n";
        errs() << space + "        if (rt < b) {\n";
        errs() << space + "            rtHistoCal(rth, b - diff, cnt);\n";
        errs() << space + "            break;\n";
        errs() << space + "        }\n";
        errs() << space + "    }\n";
        errs() << space + "}\n";
        errs() << space + "else {\n";
        errs() << space + "    rtHistoCal(rth, pow(2, msb-1), cnt);\n";
        errs() << space + "}\n";
        errs() << space + "return;\n";
        errs() << "}\n";
#ifdef REFERENCE_GROUP  
        errs() << "\n";
        errs() << "void refSubBlkRT(map<string, map<uint64_t, double>> &rth, uint64_t rt, double cnt, string ref) {\n";
        errs() << space + "int msb = 0;\n";
        errs() << space + "int tmp_rt = rt;\n";
        errs() << space + "while(tmp_rt != 0) {\n";
        errs() << space + "    tmp_rt = tmp_rt / 2;\n";
        errs() << space + "    ++msb;\n";
        errs() << space + "}\n";
        errs() << space + "if (msb >= BIN_SIZE) {\n";
        errs() << space + "    int diff = (pow(2, msb) - pow(2, msb-1)) / BIN_SIZE;\n";
        errs() << space + "    for (int b = pow(2, msb-1); b <= pow(2, msb); b+=diff) {\n";
        errs() << space + "        if (rt < b) {\n";
        errs() << space + "            refRTHistoCal(rth, b - diff, cnt, ref);\n";
        errs() << space + "            break;\n";
        errs() << space + "        }\n";
        errs() << space + "    }\n";
        errs() << space + "}\n";
        errs() << space + "else {\n";
        errs() << space + "    refRTHistoCal(rth, pow(2, msb-1), cnt, ref);\n";
        errs() << space + "}\n";
        errs() << space + "return;\n";
        errs() << "}\n";
#endif
        return;
    }
    
    
#ifdef DumpRTMR
    void AccLevelUISamplingCodeGen_ref::rtDumpGen() {
        
        errs() << "void rtDump() {\n";
        errs() << "    cout << \"Start to dump reuse time histogram\\n\";\n";
        errs() << "    for (map<uint64_t, double>::iterator it = RT.begin(), eit = RT.end(); it != eit; ++it) {\n";
        errs() << "        cout << it->first << \", \" << it->second << \"\\n\";\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
#ifdef REFERENCE_GROUP
        errs() << "\n";
        errs() << "void refRTDump() {\n";
        errs() << "    for (map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
        errs() << "        cout << \"Start to dump reuse time histogram for \" << it->first << \"\\n\";\n";
        errs() << "        for (map<uint64_t, double>::iterator iit = it->second.begin(), eiit = it->second.end(); iit != eiit; ++iit) {\n";
        errs() << "            cout << it->first << \",\" << iit->first << \",\" << iit->second << endl;\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
#endif  
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::rtToMRGen() {
        
        errs() << "void RTtoMR_AET() {\n";
        
        string space = "    ";
        errs() << space + "std::map<uint64_t, double> P;\n";
        errs() << space + "double total_num_RT = 0;\n";
        errs() << space + "uint64_t max_RT = 0;\n";
        
        errs() << space + "for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {\n";
        errs() << space + "    total_num_RT += it->second;\n";
        errs() << space + "    if (max_RT < it->first) {\n";
        errs() << space + "        max_RT = it->first;\n";
        errs() << space + "    }\n";
        errs() << space + "}\n";
        
        errs() << space + "double accumulate_num_RT = 0;\n";
        errs() << space + "for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {\n";
        errs() << space + "    P[it->first] = accumulate_num_RT / total_num_RT;\n";
        errs() << space + "    accumulate_num_RT += it->second;\n";
        errs() << space + "}\n";
        
        errs() << space + "P[0] = 1;\n";
        
        errs() << space + "double sum_P = 0;\n";
        errs() << space + "uint64_t t = 0;\n";
        errs() << space + "uint64_t prev_t = 0;\n";
//        errs() << space + "for (uint64_t c = 0; c <= max_RT; c++) {\n";
// 128M cache, 32B cache line
//        errs() << space + "for (uint64_t c = 0; c <= max_RT && c <= 1677216; c++) {\n";
// 20MB
        errs() << space + "for (uint64_t c = 0; c <= max_RT && c <= 327680; c++) {\n";
        errs() << space + "    while (sum_P < c && t <= max_RT) {\n";
        errs() << space + "        if (P.find(t) != P.end()) {\n";
        errs() << space + "            sum_P += P[t];\n";
        errs() << space + "            prev_t = t;\n";
        errs() << space + "        } else {\n";
        errs() << space + "            sum_P += P[prev_t];\n";
        errs() << space + "        }\n";
        errs() << space + "        t++;\n";
        errs() << space + "    }\n";
        
        errs() << space + "    MR[c] = P[prev_t];\n";
        errs() << space + "}\n";
        
        errs() << space + "return;\n";
        
        errs() << "}\n";
        
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::mrDumpGen() {
        
        errs() << "void dumpMR() {\n";
        
        errs() << "    cout << \"miss ratio\" << endl;\n";
        
        errs() << "    std::map<uint64_t, double>::iterator it1 = MR.begin();\n";
        errs() << "    std::map<uint64_t, double>::iterator it2 = MR.begin();\n";
        
        errs() << "    while(it1 != MR.end()) {\n";
        errs() << "        while(1) {\n";
        errs() << "            std::map<uint64_t, double>::iterator it3 = it2;\n";
        errs() << "            ++it3;\n";
        errs() << "            if (it3 == MR.end()) {\n";
        errs() << "                break;\n";
        errs() << "            }\n";
        errs() << "            if (it1->second - it3->second < 0.00001) {\n";
        errs() << "                ++it2;\n";
        errs() << "            } else {\n";
        errs() << "                break;\n";
        errs() << "            }\n";
        errs() << "        }\n";
        errs() << "        cout << it1->first << \", \" << it1->second << endl;\n";
        errs() << "        if (it1 != it2) {\n";
        errs() << "            cout << it2->first << \", \" << it2->second << endl;\n";
        errs() << "        }\n";
        errs() << "        it1 = ++it2;\n";
        errs() << "        it2 = it1;\n";
        errs() << "    }\n";
        
        errs() << "    return;\n";
        errs() << "}\n";
        
        return;
    }
#elif defined(DumpRefLease)
    void AccLevelUISamplingCodeGen_ref::accessRatioCalGen() {
        errs() << "void accessRatioCal() {\n";
        errs() << "    double total_access_cnt = 0;\n";
        errs() << "\n";
        errs() << "    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {\n";
        errs() << "        for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {\n";
        errs() << "            total_access_cnt += ri_it->second;\n";
        errs() << "        }\n";
        errs() << "    }\n";
            
        errs() << "    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {\n";
        errs() << "        double ref_access_cnt = 0;\n";
        errs() << "        for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {\n";
        errs() << "           ref_access_cnt += ri_it->second;\n";
        errs() << "       }\n";
        errs() << "        sampledCnt[ref_it->first] = ref_access_cnt;\n";
        errs() << "        accessRatio[ref_it->first] = ref_access_cnt / total_access_cnt;\n";
        errs() << "    }\n";
        errs() << "}\n";
    }
    
    void AccLevelUISamplingCodeGen_ref::initHitsCostsGen() {
        errs() << "void initHitsCosts() {\n";
            
        errs() << "    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {\n";
        errs() << "        hits[ref_it->first] = new map<uint64_t, double>;\n";
        errs() << "        costs[ref_it->first] = new map<uint64_t, double>;\n";
        errs() << "        (*hits[ref_it->first])[0] = 0;\n";
        errs() << "        uint64_t total_hits = 0;\n";
        errs() << "        (*costs[ref_it->first])[0] = 0;\n";
        errs() << "        uint64_t total_cnt = 0;\n";
        errs() << "        for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {\n";
        errs() << "            total_cnt += ri_it->second;\n";
        errs() << "        }\n";
        errs() << "        uint64_t pre_lease = 0;\n";
        errs() << "        uint64_t pre_cost = 0;\n";
        errs() << "        for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {\n";
        errs() << "            total_hits += ri_it->second;\n";
        errs() << "            (*hits[ref_it->first])[ri_it->first] = total_hits;\n";
        
        errs() << "            (*costs[ref_it->first])[ri_it->first] =  pre_cost + (ri_it->first - pre_lease) * total_cnt;\n";
        errs() << "            total_cnt -= ri_it->second;\n";
        errs() << "            pre_cost = (*costs[ref_it->first])[ri_it->first];\n";
        errs() << "            pre_lease = ri_it->first;\n";
        errs() << "        }\n";
        errs() << "    }\n";
            
        errs() << "}\n";
    }
    void AccLevelUISamplingCodeGen_ref::getPPUCGen() {
        errs() << "double getPPUC(uint64_t ref_id, uint64_t oldLease, uint64_t newLease) {\n";
            
        errs() << "    if (hits.find(ref_id) == hits.end() || costs.find(ref_id) == costs.end()) {\n";
        errs() << "        cout << \"No such ref for hits/costs\" << endl;\n";
        errs() << "        return -1;\n";
        errs() << "    }\n";
        errs() << "    if (hits[ref_id]->find(newLease) == hits[ref_id]->end() || costs[ref_id]->find(newLease) == costs[ref_id]->end()) {\n";
        errs() << "        cout << \"No RI/Newlease \" << newLease << \" for ref \" << ref_id << endl;\n";
        errs() << "        return -1;\n";
        errs() << "    }\n";
            
        errs() << "    if (hits[ref_id]->find(oldLease) == hits[ref_id]->end() || costs[ref_id]->find(oldLease) == costs[ref_id]->end()) {\n";
        errs() << "        if (hits[ref_id]->find(oldLease) == hits[ref_id]->end()) {\n";
        errs() << "            cout << \"No hits for Oldlease \" << oldLease << \" for ref \" << ref_id << endl;\n";
        errs() << "        }\n";
        errs() << "        if (costs[ref_id]->find(oldLease) == costs[ref_id]->end()) {\n";
        errs() << "            cout << \"No costs for Oldlease \" << oldLease << \" for ref \" << ref_id << endl;\n";
        errs() << "        }\n";
        errs() << "        return -1;\n";
        errs() << "    }\n";
            
        errs() << "    return double((*hits[ref_id])[newLease] - (*hits[ref_id])[oldLease]) / ((*costs[ref_id])[newLease] - (*costs[ref_id])[oldLease]);\n";
        errs() << "}\n";
    }
    void AccLevelUISamplingCodeGen_ref::getMaxPPUCGen() {
        errs() << "void getMaxPPUC(bool*finished, uint64_t* ref_to_assign, uint64_t* newLease) {\n";
            
        errs() << "    double maxPPUC = -1;\n";
        errs() << "    uint64_t bestRef = -1;\n";
        errs() << "    uint64_t bestLease = -1;\n";
            
        errs() << "    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {\n";
        errs() << "        for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {\n";
        errs() << "            if (ri_it->first > Lease[ref_it->first]) {\n";
        errs() << "                double ppuc = getPPUC(ref_it->first, Lease[ref_it->first], ri_it->first);\n";
        errs() << "                if (ppuc > maxPPUC) {\n";
        errs() << "                    maxPPUC = ppuc;\n";
        errs() << "                    bestRef = ref_it->first;\n";
        errs() << "                    bestLease = ri_it->first;\n";
        errs() << "                }\n";
        errs() << "            }\n";
        errs() << "        }\n";
        errs() << "    }\n";
            
        errs() << "    if (maxPPUC != -1) {\n";
        errs() << "        *finished = false;\n";
        errs() << "        *ref_to_assign = bestRef;\n";
        errs() << "        *newLease = bestLease;\n";
        errs() << "    } else {\n";
        errs() << "        *finished = true;\n";
        errs() << "    }\n";
            
        errs() << "    return;\n";
        errs() << "}\n";
    }
	void AccLevelUISamplingCodeGen_ref::DumpRIGen() {
        errs() << "void dumpRI() {\n";
        errs() << "    uint64_t total_number_of_ri = 0;\n";
        errs() << "    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {\n";
        errs() << "        std::set<uint64_t> riset;\n";
        errs() << "        for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {\n";
        errs() << "            cout << \"Ref \" << ref_it->first << \" RI \" << ri_it->first << \" CNT \" << ri_it->second << endl;\n";
        errs() << "            riset.insert(ri_it->first);\n";
        errs() << "        }\n";
        errs() << "        cout << \"Ref \" << ref_it->first << \" RISETSIZE \" << riset.size() << endl;\n";
        errs() << "        total_number_of_ri += riset.size();\n";
        errs() << "    }\n";
        errs() << "    cout << \"Average RISETSIZE for each reference \" << double(total_number_of_ri) / RI.size() << endl;\n";
        errs() << "}\n";
	}
    void AccLevelUISamplingCodeGen_ref::RLGen() {
        errs() << "void RL_main(uint64_t CacheSize) {\n";
        errs() << "    initHitsCosts();\n";
        errs() << "    accessRatioCal();\n";
        errs() << "    double totalCost = 0;\n";
        errs() << "    double totalHitRatio = 0;\n";
        errs() << "    double targetCost = CacheSize;\n";
        errs() << "    dumpRI();\n";
        errs() << "    while(true) {\n";
        errs() << "        bool finished = false;\n";
        errs() << "        uint64_t ref_to_assign;\n";
        errs() << "        uint64_t newLease;\n";
        errs() << "        getMaxPPUC(&finished, &ref_to_assign, &newLease);\n";
        errs() << "        if (finished == false) {\n";
        errs() << "            totalCost += ((*costs[ref_to_assign])[newLease] - (*costs[ref_to_assign])[Lease[ref_to_assign]]) / sampledCnt[ref_to_assign] * accessRatio[ref_to_assign];\n";
        errs() << "            totalHitRatio += ((*hits[ref_to_assign])[newLease] - (*hits[ref_to_assign])[Lease[ref_to_assign]]) / sampledCnt[ref_to_assign] * accessRatio[ref_to_assign];\n";
        errs() << "            Lease[ref_to_assign] = newLease;\n";
        errs() << "            cout << \"Assign lease \" << newLease << \" to ref \" << ref_to_assign << \" avg cache size \" << totalCost  << \" miss ratio \" << 1 - totalHitRatio << endl;\n";
        errs() << "        } else {\n";
        errs() << "            break;\n";
        errs() << "        }\n";
        errs() << "        if (totalCost < targetCost && targetCost != 0) {\n";
        errs() << "            break;\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
    }
#endif

#if defined(REFERENCE_GROUP)
    void AccLevelUISamplingCodeGen_ref::rtMergeGen() {
        errs() << "/* Merge the refRT to RT */\n";
        errs() << "void rtMerge() {\n";
        errs() << "    for(map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
        errs() << "        for (map<uint64_t, double>::iterator iit = it->second.begin(); iit != it->second.end(); ++iit) {\n"; 
#if defined(UNIFORM_SMOOTHING) || defined(GAUSSIAN_SMOOTHING)
        errs() << "            rtHistoCal(RT, iit->first, iit->second);\n";
#else
        errs() << "            subBlkRT(RT, iit->first, iit->second);\n";
#endif
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        return;
    }
#endif
#if defined(UNIFORM_SMOOTHING) || defined(GAUSSIAN_SMOOTHING)
    void AccLevelUISamplingCodeGen_ref::filterArrayAccesses() {
        for(map<Instruction*, vector<string>>::iterator mit = arrayAccessVariable.begin(); mit != arrayAccessVariable.end(); ++mit) {
            errs() << "/* Array " << arrayName[mit->first] << "\t";
            for(vector<string>::iterator vit = (mit->second).begin(); vit != (mit->second).end(); ++vit) {
                errs() << *vit << " ";
            }
            errs() << "*/ \n";
            bool done = true;
            for(vector<string>::iterator vit = (mit->second).begin(); vit != (mit->second).end(); ++vit) {
                for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = outloops.begin(); it != outloops.end(); ++it) {
                    // if all induction variable is outermost loop independent
                    // those references should be special considered.
                    if (indvName[(*(*it)->LIS->IDV)[0]] == *vit) {
                        done = false;
                        break;
                    }
                }
                if (!done) { 
                    break;
                }
            }
            if (done) { outMostIndependentArrayRef.push_back(mit->first); }
        }
        return;
    }
#endif
#if defined(UNIFORM_SMOOTHING)
    void AccLevelUISamplingCodeGen_ref::UniformDistrGen() {
        errs() << "/* Smoothing the per-reference RTHisto Uniformly. Equally splite the RT to a range of bins */\n";
        errs() << "void uniform_smoothing(map<uint64_t, double> &rth, bool enable, bool scale) {\n";
        errs() << "    map<uint64_t, double> tmp;\n";
        // errs() << "    for(map<uint64_t, double>::iterator it = RT.begin(); it != RT.end(); ++it) {\n";
        // errs() << "        uint64_t mu = it->first;\n";
        // errs() << "        /* Do uniform distribution for small (rt < 32) only */\n";
        // errs() << "        if (it->second <= 0.0) { continue; }\n";
        // errs() << "        uint64_t start_b = log2(mu) <= log2(THREAD_NUM) ? 0 : (log2(mu) - log2(THREAD_NUM));\n";
        // errs() << "        uint64_t end_b = (log2(mu) + log2(THREAD_NUM));\n";
        // errs() << "        double split_val = it->second / (end_b - start_b + 1);\n";
        // errs() << "        for(int b = start_b; b <= end_b; b++) {\n";
        // errs() << "            if (b >= 0) {\n";
        // errs() << "                rtHistoCal(tmp, pow(2.0, b), split_val);\n";
        // errs() << "            }\n";
        // errs() << "        }\n";
        // errs() << "    }\n";
        // errs() << "    for(map<uint64_t, double>::iterator it = tmp.begin(); it != tmp.end(); ++it) {\n";
        // errs() << "        RT[it->first] = it->second;\n";
        // errs() << "    }\n";
        // errs() << "    return;\n";
        // errs() << "}\n";
        // errs() << "    if (!enable) {\n";
        // errs() << "        for(map<uint64_t, double>::iterator it = rth.begin(); it != rth.end(); ++it) {\n";
        // errs() << "            rtHistoCal(rth, it->first, it->second);\n";
        // // errs() << "            subBlkRT(rth, it->first, it->second);\n";
        // // errs() << "        RT[it->first] = it->second;\n";
        // errs() << "        }\n";
        // errs() << "        return;\n";
        // errs() << "    }\n";
        errs() << "    double sum = accumulate(begin(rth), end(rth), 0.0, [](const double previous, const pair<uint64_t, double>& p) { return previous + p.second; });\n";
        errs() << "    double sum_P = 0.0;\n";
        errs() << "    for(map<uint64_t, double>::iterator it = rth.begin(); it != rth.end(); ++it) {\n";
        // errs() << "        if (!enable) {\n";
        // errs() << "            rtHistoCal(tmp, it->first, it->second);\n";
        // errs() << "            continue;\n";
        // errs() << "        }\n";
        errs() << "        uint64_t mu = it->first;\n";
        errs() << "        sum_P += it->second;\n";
        errs() << "        /* Do uniform distribution for all ri, distribute from ri / THREAD_NUM to ri * THREAD_NUM  */\n";
        errs() << "        if (it->second <= 0.0) { continue; }\n";
        // errs() << "        if (sum_P >= 0.9 * sum ) { \n";
        // errs() << "            rtHistoCal(tmp, mu, it->second);\n";
        // errs() << "            continue;\n";
        // errs() << "        }\n";
        errs() << "        if (mu < THREAD_NUM) { mu = mu * THREAD_NUM; }\n";
        errs() << "        uint64_t start_b = (mu / THREAD_NUM) >= 1 ? mu / THREAD_NUM : 1;\n";
        errs() << "        uint64_t end_b = mu * THREAD_NUM;\n";
        errs() << "        double split_val = it->second / (end_b - start_b + 1);\n";
        errs() << "        for (int b = start_b; b <= end_b; b++) {\n";
        errs() << "            rtHistoCal(tmp, b, split_val);\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    rth.clear();\n";
        errs() << "    for(map<uint64_t, double>::iterator it = tmp.begin(); it != tmp.end(); ++it) {\n";
        errs() << "        subBlkRT(rth, it->first, it->second);\n";
        // errs() << "        rtHistoCal(rth, it->first, it->second);\n";
        // errs() << "        RT[it->first] = it->second;\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        return;
    }
#ifdef REFERENCE_GROUP
    void AccLevelUISamplingCodeGen_ref::GroupUniformDistrGen(string space) {
        filterArrayAccesses();
        errs() << "void group_uniform_smoothing() {\n";
        double sensitive_ratio = (double)outMostIndependentArrayRef.size() / refNumber.size();
        if (sensitive_ratio >= 0.375) {
            errs() << space << "vector<string> vec = {";
            for(vector<Instruction*>::iterator it = outMostIndependentArrayRef.begin(); it != outMostIndependentArrayRef.end(); ++it) {
                errs() << "\"" << arrayName[*it] + std::to_string(refNumber[*it]) << "\"";
                if (it != outMostIndependentArrayRef.end()-1) { errs() << ", "; }
            }
            errs() << " };\n";
        } 
        errs() << space << "for (map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
        /* Cold Program */
        errs() << space << "    uniform_smoothing(it->second, true, true);";
        /*
        if (sensitive_ratio == 0.0) {
            errs() << space << "uniform_smoothing(it->second, false, false);\n";
        } else if (sensitive_ratio < 0.375) {
            errs() << space << "uniform_smoothing(it->second, true, true);\n";
        } else {
            errs() << space << "uniform_smoothing(it->second, find(vec.begin(), vec.end(), it->first) != vec.end(), true);\n";
        }
        */
        errs() << space << "}\n";
        errs() << "}\n";
    }
#endif
#elif defined(GAUSSIAN_SMOOTHING)
    void AccLevelUISamplingCodeGen_ref::GaussianDistrGen() {
        errs() << "/* Smoothing the per-reference RTHisto based on Gasussian */\n";
        errs() << "void gaussian_smoothing(map<uint64_t, double> &rth, bool scale, double sigma, bool enable) {\n";
        errs() << "    map<uint64_t, double> tmp;\n";
        errs() << "    double sum = accumulate(begin(rth), end(rth), 0.0, [](const double previous, const pair<uint64_t, double>& p) { return previous + p.second; });\n";
        errs() << "    double sum_P = 0.0;\n";
        errs() << "    for(map<uint64_t, double>::iterator it = rth.begin(); it != rth.end(); ++it) {\n";
        // errs() << "        if (!enable) {\n";
        // errs() << "            rtHistoCal(tmp, it->first, it->second);\n";
        // errs() << "            continue;\n";
        // errs() << "        }\n";
        errs() << "        sum_P += it->second;\n";
        errs() << "        uint64_t mu = it->first;\n";
        errs() << "        if (sum_P >= 0.9 * sum) {\n";
        errs() << "            rtHistoCal(tmp, mu, it->second);\n";
        errs() << "            continue;\n";
        errs() << "        }\n";
        errs() << "        if (scale && mu > THREAD_NUM) { mu = mu * THREAD_NUM; }\n";
        errs() << "        uint64_t start_b = ( mu / THREAD_NUM) >= 1 ? mu / THREAD_NUM : 1;\n";
        errs() << "        uint64_t end_b = mu;\n";
        errs() << "        if (enable) { end_b = mu * THREAD_NUM; }\n";
        errs() << "        for(uint64_t b = start_b; b <= end_b; b++) {\n";
        errs() << "            double c =  (1 /  (sqrt(2 * M_PI) * sigma));\n";
        errs() << "            double val = it->second * c * exp( -1 * pow((b - mu), 2.0) / (2 * pow(sigma, 2.0)));\n";
        errs() << "            rtHistoCal(tmp, b, val);\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    rth.clear();\n";
        errs() << "    for(map<uint64_t, double>::iterator it = tmp.begin(); it != tmp.end(); ++it) {\n";
        // errs() << "        rtHistoCal(rth, it->first, it->second);\n";
        errs() << "        subBlkRT(rth, it->first, it->second);\n";
        // errs() << "        RT[it->first] = it->second;\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        return;
    }
#ifdef REFERENCE_GROUP
    void AccLevelUISamplingCodeGen_ref::GroupGaussianDistrGen(string space) {
        filterArrayAccesses();
        errs() << "void group_gaussian_smoothing(double sigma) {\n";
        double sensitive_ratio = (double)outMostIndependentArrayRef.size() / refNumber.size();
        if (sensitive_ratio >= 0.375) {
            errs() << space << "vector<string> vec = {";
            for(vector<Instruction*>::iterator it = outMostIndependentArrayRef.begin(); it != outMostIndependentArrayRef.end(); ++it) {
                errs() << "\"" << arrayName[*it] + std::to_string(refNumber[*it]) << "\"";
                if (it != outMostIndependentArrayRef.end()-1) { errs() << ", "; }
            }
            errs() << " };\n";
        }
        errs() << space << "for (map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
        /* Cold Program */ 
        if (outMostIndependentArrayRef.size() == 0) {
            errs() << space << "/* Cold Program - No scaling */\n";
            errs() << space << space << "gaussian_smoothing(it->second, false, sigma, true);\n";
        } else if (sensitive_ratio >= 0.375) {
            errs() << space << "/* Low Sensitive Program - Scaling outmost independent array references only */\n";
            errs() << space << space << "gaussian_smoothing(it->second, find(vec.begin(), vec.end(), it->first) != vec.end(), sigma, true);\n";
        } else {
            errs() << space << "/* High Sensitive Program - Scaling all */\n";
            errs() << space << space << "gaussian_smoothing(it->second, true, sigma, true);\n";
        }
        errs() << space << "}\n";
        errs() << "}\n";
    }
#endif
#endif
    string AccLevelUISamplingCodeGen_ref::getBound(Value *bound) {
        
        if (isa<Instruction>(bound)) {
            
            Instruction *inst = cast<Instruction>(bound);
            
            switch (inst->getOpcode()) {
                case Instruction::Add:
                    return "(" + getBound(inst->getOperand(0)) + " + " + getBound(inst->getOperand(1)) + ")";
                    break;
                case Instruction::Sub:
                    return "(" + getBound(inst->getOperand(0)) + " - " + getBound(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::Mul:
                    return "(" + getBound(inst->getOperand(0)) + " * " + getBound(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::FDiv:
                case Instruction::SDiv:
                case Instruction::UDiv:
                    return "(" + getBound(inst->getOperand(0)) + " / " + getBound(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::Load:
                    return inst->getOperand(0)->getName().str();
                    break;
                case Instruction::Alloca:
                    return inst->getName().str();
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(bound)) {
            return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
        }
        return "";
    }
    
    string AccLevelUISamplingCodeGen_ref::getBound_Start(Value *bound) {
        
        if (isa<Instruction>(bound)) {
            
            Instruction *inst = cast<Instruction>(bound);
            
            switch (inst->getOpcode()) {
                case Instruction::Add:
                    return "(" + getBound_Start(inst->getOperand(0)) + " + " + getBound_Start(inst->getOperand(1)) + ")";
                    break;
                case Instruction::Sub:
                    return "(" + getBound_Start(inst->getOperand(0)) + " - " + getBound_Start(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::Mul:
                    return "(" + getBound_Start(inst->getOperand(0)) + " * " + getBound_Start(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::FDiv:
                case Instruction::SDiv:
                case Instruction::UDiv:
                    return "(" + getBound_Start(inst->getOperand(0)) + " / " + getBound_Start(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::Load:
                    return inst->getOperand(0)->getName().str() + "_Start";
                    break;
                case Instruction::Alloca:
                    return inst->getName().str() + "_Start";
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(bound)) {
            return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
        }
        return "";
    }

    string AccLevelUISamplingCodeGen_ref::getLoopInc(Value *inc) {
        if (isa<Instruction>(inc)) {
            Instruction *inst = cast<Instruction>(inc);
            switch (inst->getOpcode()) {
                case Instruction::Add:
                case Instruction::Sub:
                case Instruction::Mul:
                case Instruction::FDiv:
                case Instruction::SDiv:
                case Instruction::UDiv:
                    if (isa<ConstantInt>(inst->getOperand(0))) {
                        return getLoopInc(inst->getOperand(0));
                    } else {
                        return getLoopInc(inst->getOperand(1));
                    }
                    break;
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(inc)) {
            return to_string(dyn_cast<ConstantInt>(inc)->getValue().getSExtValue());
        }
        return "";
    }
    
    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> AccLevelUISamplingCodeGen_ref::findLoops(
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, 
        string refName, 
        int useID, 
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops,
        bool enableOPT
    ) {
        
        /* If this is a loop, add it to the list */
        if (LoopRefTree->L != NULL) {
            loops.push_back(LoopRefTree);
        }
        /* If this is a reference and it's the reference we'd like to sample, return the list */
        if (LoopRefTree->AA != NULL) {
            if (arrayName[LoopRefTree->AA] == refName) {
                if (!enableOPT && refNumber[LoopRefTree->AA] == useID) {
                    return loops;
                } else if (enableOPT && refNumber[LoopRefTree->AA] >= useID) {
                    return loops;
                }
            } else {
                return std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>();
            }
        }
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopRes;
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopTmp = findLoops(*it, refName, useID, loops, enableOPT);
                if (loopTmp.size() != 0) {
                    loopRes = loopTmp;
                }
            }
        }
        return loopRes;
    }

    uint64_t AccLevelUISamplingCodeGen_ref::getOuterMostLoopIterationSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node) {
        return ((uint64_t) stoi(getBound((*node->LIS->LB)[0].second)) - (uint64_t) stoi(getBound((*node->LIS->LB)[0].first))) / (uint64_t) stoi(getLoopInc((*node->LIS->INC)[0]));
    }

    uint64_t AccLevelUISamplingCodeGen_ref::computeIterSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {
        uint64_t iterSpace = 0;
        if (LoopRefTree->next != NULL) {
            for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                if ((*it)->L == NULL) {
                    iterSpace += 1;
                } else {
                    iterSpace += computeIterSpace(*it);
                }
            } 
            iterSpace *= getOuterMostLoopIterationSpace(LoopRefTree);
        }
        return iterSpace;
    }
    
    /* Search result reuse (Different loop) */
    /* a is coefficient vector */
    /* x is the input matrix */
    /* y is the input result */
    /* length */
    void AccLevelUISamplingCodeGen_ref::searchReuseDifferentLoopsUpdateFuncGen() {
        errs() << "void updateCoefficient(float* a, int** x, int length, int* y) {\n";
        errs() << "    double** x_tmp = new double*[length];\n";
        errs() << "    for (int i = 0; i < length; i++) {\n";
        errs() << "        x_tmp[i] = new double[length];\n";
        errs() << "        for (int j = 0; j < length; j++) {\n";
        errs() << "            x_tmp[i][j] = (double) x[i][j];\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    double* y_tmp = new double[length];\n";
        errs() << "    for (int i = 0; i < length; i++) {\n";
        errs() << "        y_tmp[i] = y[i];\n";
        errs() << "    }\n";
        errs() << "    for (int i = 0; i < length; i++) {\n";
        errs() << "        double maxEl = abs(x_tmp[i][i]);\n";
        errs() << "        int maxRow = i;\n";
        errs() << "        for (int k=i+1; k< length; k++) {\n";
        errs() << "            if (abs(x_tmp[k][i]) > maxEl) {\n";
        errs() << "                maxEl = abs(x_tmp[k][i]);\n";
        errs() << "                maxRow = k;\n";
        errs() << "            }\n";
        errs() << "        }\n";
        errs() << "        for (int k=i; k< length;k++) {\n";
        errs() << "            double tmp = x_tmp[maxRow][k];\n";
        errs() << "            x_tmp[maxRow][k] = x_tmp[i][k];\n";
        errs() << "            x_tmp[i][k] = tmp;\n";
        errs() << "        }\n";
        errs() << "        double tmp = y_tmp[maxRow];\n";
        errs() << "        y_tmp[maxRow] = y_tmp[i];\n";
        errs() << "        y_tmp[i] = tmp;\n";
        errs() << "        for (int k=i+1; k< length; k++) {\n";
        errs() << "            double c = -x_tmp[i][k]/x_tmp[i][i];\n";
        errs() << "            for (int j=i; j< length; j++) {\n";
        errs() << "                if (i==j) {\n";
        errs() << "                    x_tmp[j][k] = 0;\n";
        errs() << "                } else {\n";
        errs() << "                    x_tmp[j][k] += c * x_tmp[j][i];\n";
        errs() << "                }\n";
        errs() << "            }\n";
        errs() << "            y_tmp[k] += c * y_tmp[i];\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    for (int i=length-1; i>=0; i--) {\n";
        errs() << "        a[i] = y_tmp[i]/x_tmp[i][i];\n";
        errs() << "        for (int k=i-1; k>=0; k--) {\n";
        errs() << "            y_tmp[k] -= x_tmp[i][k] * a[i];\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "}\n";
        
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::searchReuseDifferentLoopsCalFuncGen() {
        
        errs() << "int calWithCoefficient(float* a, int *x, int length) {\n";
        errs() << "    float tmp = 0;\n";
        errs() << "    for (int i = 0; i < length; i++) {\n";
        errs() << "        tmp += a[i] * x[i];\n";
        errs() << "    }\n";
        errs() << "    tmp += a[length];\n";
        errs() << "    return (int) tmp;\n";
        errs() << "}\n";
        
        return;
    }
    
    bool AccLevelUISamplingCodeGen_ref::searchReuseDifferentLoopsInitGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag, std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> currentLoops, string space) {
        
        if (loops.size() != 0 && GenFlag == false) {
            if (LoopRefTree == loops[0]) {
                GenFlag = true;
            }
        }
        
        if (GenFlag == true) {
            if (LoopRefTree->AA != NULL) {
                if (arrayName[LoopRefTree->AA] == refName) {
                    
                    errs() << space + "/* init for diff loop check */\n";
                    int numRecord = loops.size() + 1;
                    errs() << space + "int* prev_cnt_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " = new int[" + std::to_string(numRecord) +"];\n";
                    errs() << space + "std::fill(prev_cnt_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + ", prev_cnt_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " + " + std::to_string(numRecord) + ", -1);\n";
                    
                    errs() << space + "int** x_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " = new int*["+ std::to_string(numRecord) +"];\n";
                    errs() << space + "float* a_cnt_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " = new float["+ std::to_string(numRecord) +"];\n";
                    errs() << space + "int* curr_v_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " = new int["+ std::to_string(numRecord) +"];\n";
                    int cntTmp = 0;
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        errs() << space + "int* prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " = new int[" + std::to_string(numRecord) +"];\n";
                        errs() << space + "std::fill(prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + ", prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " + " +  std::to_string(numRecord) + ", -1);\n";
                        errs() << space + "x_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "[" + std::to_string(cntTmp) + "] = prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + ";\n";
                        cntTmp++;
                    }
                    
                    errs() << space + "int * prev_one_Start_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " = new int[" + std::to_string(numRecord) +"];\n";
                    errs() << space + "std::fill(prev_one_Start_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + ", prev_one_Start_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " + " + std::to_string(numRecord) + ", 1);\n";
                    errs() << space + "x_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "[" + std::to_string(cntTmp) + "] = prev_one_Start_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + ";\n";
                    
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        errs() << space + "int prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "[" + std::to_string(numRecord) +"];\n";
                        errs() << space + "float* a_" + indvName[(*(*it)->LIS->IDV)[0]] + "_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " = new float["+ std::to_string(numRecord) +"];\n";
                        errs() << space + "std::fill(prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + ", prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " + " +  std::to_string(numRecord) + ", -1);\n";
                    }
                    errs() << space + "int prev_cnt_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "_idx = 0;\n";
                }
            }
        }
        
        if (LoopRefTree->next != NULL) {
            if (LoopRefTree == loops.back()) {
                return GenFlag;
            }
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops_New = currentLoops;
                if (LoopRefTree->L != NULL) {
                    currentLoops_New.push_back(LoopRefTree);
                }
                GenFlag = searchReuseDifferentLoopsInitGen(*it, GenFlag, refName, useID, loops, currentLoops_New, space);
            }
        }
        
        return GenFlag;
    }
    
    bool AccLevelUISamplingCodeGen_ref::searchReuseDifferentLoopsBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag, std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> currentLoops, string space) {
        
        if (loops.size() != 0 && GenFlag == false) {
            if (LoopRefTree == loops[0]) {
                GenFlag = true;
            }
        }
        
        if (GenFlag == true) {
            if (LoopRefTree->AA != NULL) {
                if (arrayName[LoopRefTree->AA] == refName) {
                 
                    /* genrate if statement to check coefficients are initalized or not */
                    errs() << space + "if ( ";
                    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator firstIndv = loops.begin();
                    for (unsigned long i = 0; i < currentLoops.size(); i++) {
                        errs() << "prev_" + indvName[(*(*firstIndv)->LIS->IDV)[0]] + "_Start_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "[" + std::to_string(i) + "] != -1";
                        if (i + 1 != currentLoops.size()) {
                            errs() << " && ";
                        }
                    }
                    errs() << ") {\n";

                    /* generate assignment statements to init curr_v_ */
                    for (unsigned long i = 0; i < loops.size(); i++) {
                        errs() << "curr_v_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "[" + std::to_string(i) + "] = " + indvName[(*(loops[i])->LIS->IDV)[0]] + "_Start;\n";
                    }
                    errs() << "curr_v_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "[" + std::to_string(loops.size()) + "] = 1;\n";
                    
                    /* generate if statement to check prediction access the same location or not */
                    errs() << space + "    if ( calAddr" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "( ";
                    for (unsigned long i = 0; i < currentLoops.size(); i++) {
                        errs() << "calWithCoefficient( ";
                        for (unsigned long j = 0; j < loops.size(); j++) {
                            
                        }
                        errs() << "a_" + indvName[(*(currentLoops[i])->LIS->IDV)[0]] + "_" + refName + std::to_string(refNumber[LoopRefTree->AA]);
                        errs() << " , ";
                        errs() << "curr_v_" + refName + std::to_string(refNumber[LoopRefTree->AA]);
                        errs() << " , ";
                        errs() << std::to_string(loops.size());
                        errs() << ")";
                        if (i + 1 != currentLoops.size()) {
                            errs() << ",";
                        }
                    }
                    errs() << + ") == ";

                    errs() << "calAddr" + refName + std::to_string(useID) + "( ";
                    for (unsigned long i = 0; i < loops.size(); i++) {
                        errs() << indvName[(*(loops[i])->LIS->IDV)[0]] + "_Start";
                        if (i + 1 != loops.size()) {
                            errs() << " , ";
                        }
                    }
                    errs() << ")";
                    errs() << " ) {\n";

#ifdef PROFILE_SEARCH_REUSE
                    errs() << space + "        pred_num++;\n";
                    errs() << space + "        pred_dl_num++;\n";
                    errs() << space + "        pred = ";
                    errs() << " calWithCoefficient(";
                    errs() << "a_cnt_" + refName + std::to_string(refNumber[LoopRefTree->AA]);
                    errs() << " , ";
                    errs() << "curr_v_" + refName + std::to_string(refNumber[LoopRefTree->AA]);
                    errs() << " , ";
                    errs() << std::to_string(loops.size());
                    errs() << ");\n";
#endif
                    
                    /* predict cnt and accumulate */
                    errs() << space + "        rtHistoCal(";
                    errs() << " calWithCoefficient(";
                    errs() << "a_cnt_" + refName + std::to_string(refNumber[LoopRefTree->AA]);
                    errs() << " , ";
                    errs() << "curr_v_" + refName + std::to_string(refNumber[LoopRefTree->AA]);
                    errs() << " , ";
                    errs() << std::to_string(loops.size());
                    errs() << ")";
#ifdef DumpRefLease
                    errs() << ", " + std::to_string(refNumber[LoopRefTree->AA]);
#endif
                    errs() << " );\n";

#ifdef PROFILE_SEARCH_REUSE
                    errs() << space + "        goto PROFILE;\n";
#else
                    errs() << space + "        goto EndSample;\n";
#endif
                    errs() << space + "    }\n";
                    errs() << space + "}\n";

                }
            }
        }
        
        if (LoopRefTree->next != NULL) {
            if (LoopRefTree == loops.back()) {
                return GenFlag;
            }
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops_New = currentLoops;
                if (LoopRefTree->L != NULL) {
                    currentLoops_New.push_back(LoopRefTree);
                }
                GenFlag = searchReuseDifferentLoopsBodyGen(*it, GenFlag, refName, useID, loops, currentLoops_New, space);
            }
        }
        
        return GenFlag;
    }
    

    /* Search result reuse (Same loop) */
    void AccLevelUISamplingCodeGen_ref::searchReuseSameLoopInitGen(std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, string space) {
        
        //errs() << "Search result reuse: (ref name " << refName << " ) ( ID " << useID << " ) ( numOfLoops " << loops.size() << " )\n";

        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator ait = loops.back()->next->begin(), eit = loops.back()->next->end(); ait != eit; ++ait) {
            if ((*ait)->isThreadNode) { continue; }
            if ((*ait)->AA != NULL) {
                if (arrayName[(*ait)->AA] == refName) {
                    errs() << space + "uint64_t prev_cnt_" + refName + std::to_string(refNumber[(*ait)->AA]) + " = -1;\n";
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        errs() << space + "uint64_t prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName + std::to_string(refNumber[(*ait)->AA]) + " = -1;\n";
                        errs() << space + "uint64_t prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName + std::to_string(refNumber[(*ait)->AA]) + " = -1;\n";
                    }
                }
            }
        }
    
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::searchReuseSameLoopBodyGen(std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops,  string space) {
        
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator ait = loops.back()->next->begin(), eit = loops.back()->next->end(); ait != eit; ++ait) {
            if ((*ait)->isThreadNode) { continue; }
            if ((*ait)->AA != NULL) {
                if (arrayName[(*ait)->AA] == refName) {
                    errs() << space + "if ( prev_cnt_" + refName + std::to_string(refNumber[(*ait)->AA]) + " != -1) {\n";
                    errs() << space + "    if ( calAddr" + refName + std::to_string(refNumber[(*ait)->AA]) + "( ";
                    string tmp = "";
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            
                            tmp += indvName[(*(*it)->LIS->IDV)[i]] + "_Start";
                            tmp += " - prev_" + indvName[(*(*it)->LIS->IDV)[i]] + "_Start_" + refName + std::to_string(refNumber[(*ait)->AA]);
                            tmp += " + prev_" + indvName[(*(*it)->LIS->IDV)[i]] + "_End_" + refName + std::to_string(refNumber[(*ait)->AA]);
                            tmp += ", ";
                        }
                    }
                    if (loops.size() != 0) {
                        tmp.pop_back();
                        tmp.pop_back();
                    }
                    
                    errs() << tmp + ")";
                    errs() << " == ";
                    
                    errs() << "calAddr" + refName + std::to_string(useID);
                    errs() << "(";
                    
                    tmp = "";
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            tmp += indvName[(*(*it)->LIS->IDV)[i]];
                            tmp += "_Start";
                            tmp += ", ";
                        }
                    }
                    if (loops.size() != 0) {
                        tmp.pop_back();
                        tmp.pop_back();
                    }
                    errs() << tmp + ")) {\n";

#ifdef PROFILE_SEARCH_REUSE
                    errs() << space + "        pred_num++;\n";
                    errs() << space + "        pred_sl_num++;\n";
                    errs() << space + "        pred = " + "prev_cnt_" + refName + std::to_string(refNumber[(*ait)->AA]) + ";\n";
#endif
                    
                    errs() << space + "        rtHistoCal(prev_cnt_" + refName + std::to_string(refNumber[(*ait)->AA]);
#ifdef DumpRefLease
                    errs() << ", " + std::to_string(refNumber[(*ait)->AA]);
#endif
                    errs() << ");\n";


#ifdef PROFILE_SEARCH_REUSE
                    errs() << space + "        goto PROFILE;\n";
#else
                    errs() << space + "        goto EndSample;\n";
#endif
                    errs() << space + "    }\n";
                    errs() << space + "}\n";
                }
            }
        }
        
        return;
    }
    
    /* Generating Reuse Search */
    bool AccLevelUISamplingCodeGen_ref::refRTSearchGen(
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, 
        bool GenFlag,
        bool isFirstNestLoop, 
        std::string refName, 
        int useID, 
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, 
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, 
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> sampleIDVs, 
        string space) {
#ifdef PSCODEGEN_DEBUG
        errs() << space + " /* " << loops.size() << " */\n";
#endif
        
        if (loops.size() != 0 && GenFlag == false) {
            if (LoopRefTree == loops[0]) {
#ifdef PSCODEGEN_DEBUG
                errs() << space + " /* GenFlag = true */\n";
#endif
                GenFlag = true;
            }
        }



#ifdef PARALLEL
        bool hasAA = false;
#endif
        if (GenFlag == true) {
#ifdef PARALLEL
            // This is the LoopRefNode
            if (!LoopRefTree->isThreadNode) {
                // check if is the out-most loops
                if (find(outloops.begin(), outloops.end(), LoopRefTree) != outloops.end()) {
                    // is out-most loops
                    // compute the per-thread chunk
                    errs() << "#ifdef DEBUG\n";
                    errs() << space + "cout << \"Count: \" << cnt << endl;\n";
                    errs() << "#endif\n";
                    errs() << space + "/* Compute the chunk size. */\n";
                    errs() << "#ifdef CHUNK_SIZE\n";
                    errs() << space + "chunk_size = CHUNK_SIZE;\n";
                    errs() << space + "chunk_num = ";
                    if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULE) {
                        errs() << "(" + getBound_Start((*LoopRefTree->LIS->LB)[0].second) + " - " + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + ") % (THREAD_NUM * chunk_size) == 0 ? " + "(" + getBound_Start((*LoopRefTree->LIS->LB)[0].second) + " - " + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + ") / (THREAD_NUM * chunk_size) : " + "(" + getBound_Start((*LoopRefTree->LIS->LB)[0].second) + " - " + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + ") / (THREAD_NUM * chunk_size) + 1;\n";
                    } else if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGT) {
                        errs() << "(" + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + " - " + getBound_Start((*LoopRefTree->LIS->LB)[0].second) + ") % (THREAD_NUM * chunk_size) == 0 ? " + "(" + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + " - " + getBound_Start((*LoopRefTree->LIS->LB)[0].second) + ") / (THREAD_NUM * chunk_size) : " + "(" + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + " - " + getBound_Start((*LoopRefTree->LIS->LB)[0].second) + ") / (THREAD_NUM * chunk_size) + 1;\n";
                    } else {
                        errs() << "\n Error in computing thread local iteration space \n";
                    }
                    errs() << "#else\n";
                    errs() << space + "chunk_num = 1;\n";
                    errs() << space + "chunk_size = ";
                    if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULE) {
                        errs() << "(" + getBound_Start((*LoopRefTree->LIS->LB)[0].second) + " - " + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + ") / THREAD_NUM;\n";
                    } else if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGT) {
                        errs() << "(" + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + " - " + getBound_Start((*LoopRefTree->LIS->LB)[0].second) + ") / THREAD_NUM;\n";
                    } else {
                        errs() << "\n Error in computing thread local iteration space \n";
                    }
                    errs() << "#endif\n";

                    /* Generate the code to compute the start chunk of the sampler */
                    string loopNum = std::to_string(loopNumber[LoopRefTree->L]);
                    errs() << space + "/* Compute the number of chunks */\n";
                    if (isFirstNestLoop) {
                        errs() << space + "c_Start = ";
                        errs() << "(" + indvName[(*LoopRefTree->LIS->IDV)[0]] + "_Start - " + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + ") / (THREAD_NUM * chunk_size);\n";
                    } else {
                        errs() << space + "c_Start = 0;\n";
                    }
                    errs() << "#ifdef DEBUG\n";
                    errs() << space + "cout << \"c_Start = \" << c_Start << \", chunk_num = \" << chunk_num << endl;\n";
                    errs() << "#endif\n";
                    /* Outer-most loops will generate chunk iterations */
                    errs() << space + "/* Generating thread local iteration space mapping code */\n";
                    errs() << space + "for (int cid = c_Start; cid < chunk_num; cid++) {\n";
                    errs() << space + "    /* Computes bound express for each thread */\n";
                    errs() << space + "    for (int t = 0; t < THREAD_NUM; ++t) {\n";
                    errs() << space + "        BLIST[t][0] =  " + getBound((*LoopRefTree->LIS->LB)[0].first) + "+ (cid * THREAD_NUM + t) * chunk_size;\n";
                    errs() << space + "        BLIST[t][1] = min(" + getBound((*LoopRefTree->LIS->LB)[0].first) + " + (cid * THREAD_NUM + t + 1) * chunk_size, " + getBound((*LoopRefTree->LIS->LB)[0].second) + ") - 1;\n";
                    errs() << "#ifdef DEBUG\n";
                    errs() << space + "        cout << \"[Thread \" << t << \"], \" << \"(\" << BLIST[t][0] << \", \"<< BLIST[t][1] << \")\" << endl;\n";
                    errs() << "#endif\n"; 
                    errs() << space + "    }\n";
                    errs() << space + "    /* Iterate within a chunk */\n";
                    if (isFirstNestLoop) {
                        errs() << space + "    ci_Start = 0;\n";
                        errs() << space + "    if (cid == c_Start) {\n";
                        errs() << space + "        ci_Start = (" + indvName[(*LoopRefTree->LIS->IDV)[0]] + "_Start - " + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + ") \% chunk_size;\n";
                        errs() << space + "    }\n";
                        errs() << space + "    int " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "_Start;\n";
                    } else {
                        errs() << space + "    ci_Start = 0;\n";
                        errs() << space + "    int " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = " + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + ";\n";
                    }
                    errs() << space + "    for ( int ci = ci_Start; ci < chunk_size; ci++) {\n";
                    errs() << space + "        if ( cid != c_Start || ci != ci_Start ) {\n";
                    errs() << space + "            " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = cid * (THREAD_NUM * chunk_size) + ci;\n";
                    errs() << space + "        }\n";
                }
#endif
                /* Generate loop */
                /* This is a Loop Node */
                if (LoopRefTree->L != NULL) {
                    string loopNum = std::to_string(loopNumber[LoopRefTree->L]);
                    /* If this is the sampled loop, the following if-cond will be satisfied and the start points of inner loops will be computed */
                    if (find(outloops.begin(), outloops.end(), LoopRefTree) == outloops.end()) {
                        errs() << space + "    int " + indvName[(*LoopRefTree->LIS->IDV)[0]] +  "LB" + loopNum + " = " + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + ";\n";
                    }
                    if (std::find(loops.begin(), loops.end(), LoopRefTree) != loops.end()) {
                        if (isFirstNestLoop && LoopRefTree != loops[0]) {
                            errs() << space + "    if ( ";
                            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                                if (LoopRefTree == *it) {
                                    break;
                                }
                                if (it != loops.begin()) {
                                    errs() << " && ";
                                }
                                errs() << indvName[(*(*it)->LIS->IDV)[0]];
                                if (find(outloops.begin(), outloops.end(), (*it)) != outloops.end()) {
                                    errs() << "LB" + to_string(loopNumber[(*it)->L]);
                                }
                                errs() << " == " + indvName[(*(*it)->LIS->IDV)[0]] + "_Start";
                            }
                            errs() << " ) {\n";
                            errs() << space + "        " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = ";
                            errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] + "_Start;\n";
                            errs() << space + "    }\n";
                        }
                    }
                    /* Parallel Interleaving*/
#ifdef PARALLEL
                    /* Generate the code to compute the iteration space for each thread, the bound should be an expression with respect to the out-most loop induction variable and thread id */
                    if (find(outloops.begin(), outloops.end(), LoopRefTree) == outloops.end()) {
                        errs() << space + "    for ( int " + indvName[(*LoopRefTree->LIS->IDV)[0]];
                        errs() << " = " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum;
                        errs() << "; ";
                        errs() << indvName[(*LoopRefTree->LIS->IDV)[0]];
                        if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULE) {
                            errs() << " <= ";
                        } else if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGE) {
                            errs() << " >= ";
                        } else if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULT) {
                            errs() << " < ";
                        } else if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGT) {
                            errs() << " > ";
                        } else {
                            errs() << "\n Error recognizing predicates \n";
                        }
                        errs() << getBound_Start((*LoopRefTree->LIS->LB)[0].second);
                        errs() << "; ";
                        /* need to take stride into consideration */
                        if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULT) {
                            
                            errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] << "=" << getBound((*LoopRefTree->LIS->INC)[0]);
                            
                        } else if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGT) {
                            
                            errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] << "=" << getBound((*LoopRefTree->LIS->INC)[0]);
                            
                        } else {
                            errs() << "\n Error in geting stride \n";
                        }
                        errs() << ") {\n";
                    }
                    if (LoopRefTree->next != NULL) {
                        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                            // contains a thread node is the node that contains a array access
                            hasAA = hasAA || (*it)->isThreadNode;
                        }
                    }
#ifdef PSCODEGEN_DEBUG
                            errs() << space + "/* Loop: " << LoopRefTree->L->getName() << "hasAA: " << hasAA << " */\n"; 
#endif

                    /* Reach the inner-most loop */
                    if (hasAA) {
                        /* Check if this is the first compare in if-cond */
                        errs() << space + "        int ";
                        /* if the reference is enclosed in the out-most loop, currentLoops is empty */
                        if (currentLoops.size() == 0) {
                            errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] + " = cid * (THREAD_NUM * chunk_size) + ci + " + getBound((*LoopRefTree->LIS->LB)[0].first) + ";\n";
                            errs() << space + "        if(" + indvName[(*LoopRefTree->LIS->IDV)[0]] + " > BLIST[0][1]) { goto EndSample; }\n";
                        } else {
                            errs() << indvName[(*(*currentLoops.begin())->LIS->IDV)[0]] + " = cid * (THREAD_NUM * chunk_size) + ci + " + getBound((*currentLoops.front()->LIS->LB)[0].first) + ";\n";
                            errs() << space + "        if(" + indvName[(*(*currentLoops.begin())->LIS->IDV)[0]] + " > BLIST[0][1]) { goto EndSample; }\n";
                        }
                        // for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        //     if (std::find(sampleIDVs.begin(), sampleIDVs.end(), (*it)) != sampleIDVs.end()) {
                        //         if (it == currentLoops.begin()) {
                        //             errs() << space + "    int " + indvName[(*(*it)->LIS->IDV)[0]] + " = cid * (THREAD_NUM * chunk_size) + ci + " + getBound((*currentLoops.front()->LIS->LB)[0].first) + ";\n";
                        //             errs() << space + "    if(" + indvName[(*(*it)->LIS->IDV)[0]] + " > BLIST[0][1]) { goto EndSample; }\n";
                        //         }
                        //         else {
                        //             errs() << " && ";
                        //         }
                        //         if (find(outloops.begin(), outloops.end(), *it) == outloops.end()) {
                        //             errs() << indvName[(*(*it)->LIS->IDV)[0]] + " != " + indvName[(*(*it)->LIS->IDV)[0]] + "_Start ";
                        //             hasPrevComp = true;
                        //         }
                        //     }
                        //     if ((*it) == (*currentLoops.rbegin()) && std::find(sampleIDVs.begin(), sampleIDVs.end(), LoopRefTree) == sampleIDVs.end()) {
                        //         errs() << " ) {\n";
                        //     }
                        // }
                        // if (std::find(sampleIDVs.begin(), sampleIDVs.end(), LoopRefTree) != sampleIDVs.end()) {
                        //     if (hasPrevComp) {
                        //         errs() << "&& ";
                        //     }
                        //     else {
                        //         errs() << space + "    if ( ";
                        //     }
                        //     errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] + " != " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "_Start) {\n";
                        // }
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "        cout << \"Iterate (\"";
                        string tmp = "";
                        for(vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(); it != currentLoops.end(); ++it) {
                            tmp += (" << " + indvName[(*(*it)->LIS->IDV)[0]] + " << \", \"");
                        }
                        if (currentLoops.size() == 0) {
                            tmp += " << " + indvName[(*LoopRefTree->LIS->IDV)[0]] + " << ";
                        }
                        errs() << tmp + "\")\" << endl;\n";
                        errs() << "#endif\n";
                        /* Generate the interleaving iteation space */
                        errs() << space + "        vector<int> v = { ";
                        tmp = "";
                        for(vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(); it != currentLoops.end(); ++it) {
                            tmp += (indvName[(*(*it)->LIS->IDV)[0]] + ", ");
                        }
                        tmp.pop_back();
                        tmp.pop_back();
                        // tmp += indvName[(*LoopRefTree->LIS->IDV)[0]];
                        errs() << tmp + " };\n";
                        errs() << space + "        /* Interleaving */\n";
                        if (currentLoops.size() == 0) {
                            errs() << space + "        t_Start = ((" + indvName[(*(sampleIDVs.front())->LIS->IDV)[0]] + " - " + getBound_Start((*LoopRefTree->LIS->LB)[0].first) + ") / chunk_size) % THREAD_NUM;\n";
                        } else {
                            errs() << space + "        t_Start = ((" + indvName[(*(sampleIDVs.front())->LIS->IDV)[0]] + " - " + getBound_Start((*currentLoops.front()->LIS->LB)[0].first) + ") / chunk_size) % THREAD_NUM;\n";
                        }
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "        cout << \"Generate interleaved iteration for (\";\n";
                        errs() << space + "        for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {\n";
                        errs() << space + "            cout << *it;\n";
                        errs() << space + "            if (it != v.end()) { cout << \", \"; }\n";
                        errs() << space + "        }\n";
                        errs() << space + "        cout << \")\" << endl;\n";
                        errs() << "#endif\n";
                        errs() << space + "        for ( int tid = t_Start; tid < THREAD_NUM; tid++) {\n"; 
                        errs() << space + "            nv[tid] = v;\n";
                        errs() << space + "            nv[tid][0] = v[0] + chunk_size * (tid - t_Start);\n";
                        // errs() << space + "            vector<int> tmp;\n";
                        // errs() << space + "            for (int vi = 0; vi < v.size(); vi++ ) {\n";
                        // errs() << space + "                if (vi == 0) {\n";
                        // errs() << space + "                    tmp.push_back(v[0] + chunk_size * (tid - t_Start));\n";
                        // errs() << space + "                } else {\n";
                        // errs() << space + "                    tmp.push_back(v[vi]);\n";
                        // errs() << space + "                }\n";
                        // errs() << space + "            }\n";
                        // errs() << space + "            if (tmp.size() > 0) { nv[tid] = tmp; }\n";
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "            cout << \"(\";\n";
                        errs() << space + "            for (vector<int>::iterator it = nv[tid].begin(); it != nv[tid].end(); it++) {\n";
                        errs() << space + "                cout << *it << \", \";\n"; 
                        errs() << space + "            }\n";
                        errs() << space + "            cout << \")\" << endl;\n";
                        errs() << "#endif\n";
                        errs() << space + "        }\n";
                    }
                }
#endif
                /* Generate addr checking instruction */
                if (LoopRefTree->AA != NULL) {
                    if (arrayName[LoopRefTree->AA] == refName) {
                        errs() << space + "    /* Remove those invalid interleaving */\n";
                        errs() << space + "    if (nv[nvi].size() <= 0) { continue; }\n";
                        errs() << space + "    if (nv[nvi][0] > BLIST[nvi][1]) { break; }\n";
                        errs() << space + "    if (cntStart == true) {\n";
                        errs() << space + "        cnt++;\n";
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "        cout  << \"[" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "]\" << ";
                        for (unsigned i = 0; i < currentLoops.size(); i++) {
                            errs() << "nv[nvi][" + std::to_string(i) + "] << ";
                            if (i != currentLoops.size() - 1) {
                                errs() << "\", \" << ";
                            }
                        }
                        errs() << "\", cnt: \" << cnt << \")\t\";\n";
                        errs() << "#endif\n";                  
                        errs() << space + "        if ( calAddr" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "( ";
                        string tmp = "";
                        for (unsigned i = 0; i < currentLoops.size(); i++) {
                            tmp += "nv[nvi][" + std::to_string(i) + "], ";
                        }
                        if (currentLoops.size() != 0) {
                            tmp.pop_back();
                            tmp.pop_back();
                        }

                        errs() << tmp + ")";
                        errs() << " == ";
                        
                        errs() << "calAddr" + refName + std::to_string(useID);
                        errs() << "(";
                        
                        tmp = "";
                        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                            for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                                tmp += indvName[(*(*it)->LIS->IDV)[i]];
                                tmp += "_Start";
                                tmp += ", ";
                            }
                        }

                        if (loops.size() != 0) {
                            tmp.pop_back();
                            tmp.pop_back();
                        }
                        errs() << tmp + ")) {\n";

                        errs() << "#if defined(INTERLEAVE_DEBUG) || defined(DEBUG)\n";
                        errs() << space + "            cout << cnt << \" " << refName << std::to_string(useID) << " (\" << ";
                        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                            for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                                errs() << indvName[(*(*it)->LIS->IDV)[i]] << "_Start";
                                if (it != loops.end()-1) {
                                    errs() << " << \",\" << ";
                                }
                            }
                        }
                        errs() << "<< \") -> " << refName << std::to_string(refNumber[LoopRefTree->AA]) << " (\" << " ;
                        tmp = "";
                        for (unsigned i = 0; i < currentLoops.size(); i++) {
                            tmp += "nv[nvi][" + std::to_string(i) + "]";
                            if (i != currentLoops.size() - 1) {
                                tmp += " << \",\" << ";
                            }
                        }
                        errs() << tmp << " << \")\" << endl;\n";


                        // errs() << space + "            cout << \"[" << refName << std::to_string(useID) << " --> " << refName << std::to_string(refNumber[LoopRefTree->AA]) << "] @ (\" << ";
                        // for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        //     for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                        //         errs() << indvName[(*(*it)->LIS->IDV)[i]] << "_Start";
                        //         if (it != loops.end()-1) {
                        //             errs() << " << \", \" << ";
                        //         }
                        //     }
                        // }
                        // errs() << " << \") --> (\" << ";
                        // tmp = "";
                        // for (unsigned i = 0; i < currentLoops.size(); i++) {
                        //     tmp += "nv[nvi][" + std::to_string(i) + "]";
                        //     if (i != currentLoops.size() - 1) {
                        //         tmp += " << \", \" << ";
                        //     }
                        // }
                        // errs() << tmp + " << \"), \" << cnt << \") \" << endl;\n";
                        // errs() << space + "            if (cnt == 2907  || cnt == 2955) { exit(1); }\n";
                        errs() << "#endif\n";

#ifdef REFERENCE_GROUP
                        // errs() << space + "            refSubBlkRT(refRT, cnt, 1.0, \"" + refName + std::to_string(useID) + "\");\n";
                        errs() << space + "            refRTHistoCal(refRT, cnt, 1.0, \"" + refName + std::to_string(useID) + "\");\n";
#else
                        // errs() << space + "            subBlkRT(RT, cnt, 1.0);\n";
                        errs() << space + "            rtHistoCal(RT, cnt, 1.0);\n";
#endif
                        errs() << space + "            goto EndSample;\n";
                        
                        errs() << space + "        }\n";
                        errs() << space + "    }\n";
                        if (useID == refNumber[LoopRefTree->AA]) {
                            errs() << space + "    if (";
                            for (unsigned i = 0; i < currentLoops.size(); i++) {
                                if (i != 0) {
                                    errs() << " && ";
                                }
                                errs() << "nv[nvi][" + std::to_string(i) + "] == " + indvName[(*currentLoops[i]->LIS->IDV)[0]] + "_Start";
                            }
                            errs() << ") { cntStart = true; }\n";
                        }
                        errs() << space + "    }\n";
#ifdef PARALLEL
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "    cout << endl;\n";
                        errs() << space + "    /* useID: " + to_string(useID) + " refNumber[LoopRefTree->AA]: " + to_string(refNumber[LoopRefTree->AA]) + " */\n";
                        errs() << "#endif\n";
#endif
                    } else {
                        errs() << space + "    if (nv[nvi].size() <= 0) { continue; }\n";
                        errs() << space + "    if (nv[nvi][0] > BLIST[nvi][1]) { break; }\n";
                        errs() << space + "    if (cntStart == true) {\n";
                        errs() << space + "        cnt++;\n";
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "        cout  << \"[" + arrayName[LoopRefTree->AA] + std::to_string(refNumber[LoopRefTree->AA]) + "]\" << ";
                        for (unsigned i = 0; i < currentLoops.size(); i++) {
                            errs() << "nv[nvi][" + std::to_string(i) + "] << ";
                            if (i != currentLoops.size() - 1) {
                                errs() << "\", \" << ";
                            }
                        }
                        errs() << "\", cnt: \" << cnt << \")\t\";\n";
                        errs() << "#endif\n";     
                        errs() << space + "    }\n";
#ifdef PARALLEL
                        errs() << space + "} // end of interleaving loop\n";
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "cout << endl;\n";
                        errs() << "#endif\n";
#endif
                    }
                }   
            } else {
                // this is a thread node, generate the interleaving loop
#ifdef PARALLEL
                errs() << space + "    /* iterate thread local iteration space mapping code after interleaving */\n";
                errs() << space + "    for (int nvi = 0; nvi < nv.size(); nvi++) {\n";
#endif
            }
        }
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops_New = currentLoops;
                if ((*it)->L != NULL) {
                    currentLoops_New.push_back(*it);
                }
                /* isFirstNestLoop will be updated when we meet a out-most loops */
                if (find(outloops.begin(), outloops.end(), *it) != outloops.end()) {
                    isFirstNestLoop = ((*it) == loops[0]);
                    /* Search the new nested loop to check if it contains references that we are currently sampling */
                    if (!isFirstNestLoop && findLoops(*it, refName, useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(), true).size() <= 0) {
                        errs() << space << "    if (cntStart == true) {\n";
                        errs() << space << "        cnt += " << computeIterSpace(*it) << ";\n";
                        errs() << space << "    }\n";
                        continue;
                    }
                }
                GenFlag = refRTSearchGen(*it, GenFlag, isFirstNestLoop, refName, useID, loops, currentLoops_New, sampleIDVs, space + "    ");
            }
            if (LoopRefTree->L != NULL && GenFlag == true) {
                if (find(outloops.begin(), outloops.end(), LoopRefTree) != outloops.end()) {
                    errs() << space + "    } // end of outer for - ci loops\n";
                    errs() << space + "} // end of outer for - cid loops\n";
                } else {
                    errs() << space + "} // end of inner for loops\n";
                }
// #ifdef PARALLEL
//                 if (find(outloops.begin(), outloops.end(), LoopRefTree) == outloops.end()) {
//                     errs() << space + "if (cntStart == true) {\n";
//                     errs() << space + "    threadLB = 0;\n";
//                     errs() << space + "}\n";
//                 }
// #endif
            }
        }
    
        return GenFlag;
    }


    void AccLevelUISamplingCodeGen_ref::refRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, string refName, int useID) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops = findLoops(LoopRefTree, refName, useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(), false);
        
        for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(); lit != loops.end(); lit++) {
            errs() << "/* for (" << indvName[(*(*lit)->LIS->IDV)[0]] << ", " << getBound_Start((*(*lit)->LIS->LB)[0].first) << ", " << getBound_Start((*(*lit)->LIS->LB)[0].second)  << ") */\n";
        }


        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> sampleIDVs;
        
#ifdef PROFILE_SEARCH_REUSE
        errs() << "    /* Generating profile counter */\n";
        errs() << "    uint64_t pred_num = 0;\n";
        errs() << "    uint64_t pred_sl_num = 0;\n";
        errs() << "    uint64_t pred_dl_num = 0;\n";
        errs() << "    uint64_t sample_num = 0;\n";
        errs() << "    uint64_t reuse_in_same_loop = 0;\n";
        errs() << "    uint64_t reuse_in_diff_loops = 0;\n";
        errs() << "    uint64_t no_reuse_num = 0;\n";
        errs() << "    uint64_t mis_pred_num = 0;\n";
        errs() << "    uint64_t pred = 0;\n";
        errs() << "    uint64_t actural = 0;\n";
#endif
        
#if SAMPLING ==2
        errs() << "    /* Generating sampling loop */\n";
        string space = "    ";
        errs() << space + "set<string> record;\n";
        
        errs() << space + "for ( int s = 0; s < ";
		if (loops.size() == 0) {
			errs() << "1";
		} else {
        	if (sampleNum.find(loops.back()) != sampleNum.end()) {
            	errs() << std::to_string(sampleNum[loops.back()]);
        	} else {
            	errs() << "ERROR in finding bounds\n";
        	}
		}
        errs() << ";) {\n";
        
        if (loops.size() == 0) {
            errs() << space + "/* Generating reuse search code */\n";
            errs() << "\n";
            refRTSearchGen(LoopRefTree, false, true, refName, useID, loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(), sampleIDVs, "    ");
            errs() << space + "s++;\n";
            errs() << space + "}\n";
            return;
        }

        errs() << "SAMPLE:\n";
        
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {

                if (lit != loops.begin()) {
                    if ((*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SGE || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SLE || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_ULE || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_UGE) {
                        errs() << space + "    if ( (" + (getBound_Start((*(*lit)->LIS->LB)[i].second) + " - " + getBound_Start((*(*lit)->LIS->LB)[i].first)) + " + 1) == 0) goto SAMPLE;\n";
                    } else if ((*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SGT || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SLT || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_ULT || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_UGT) {
                        errs() << space + "    if ( (" + (getBound_Start((*(*lit)->LIS->LB)[i].second) + " - " + getBound_Start((*(*lit)->LIS->LB)[i].first)) + ") == 0) goto SAMPLE;\n";
                    } else {
                        errs() << "\n Error in generating random sample \n";
                    }
                }
               
				/* generate sampling statement */ 
                errs() << space + "    int " + indvName[(*(*lit)->LIS->IDV)[i]] + "_Start" + " = ";
                sampleIDVs.push_back((*lit));

                if ( (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SLE || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_ULE ) {
                    errs() << "rand() % (" + (getBound_Start((*(*lit)->LIS->LB)[i].second) + " - " + getBound_Start((*(*lit)->LIS->LB)[i].first)) + " + 1) + " + getBound_Start((*(*lit)->LIS->LB)[i].first);
                } else if ((*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SLT || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_ULT) {
                    errs() << "rand() % (" + (getBound_Start((*(*lit)->LIS->LB)[i].second) + " - " + getBound_Start((*(*lit)->LIS->LB)[i].first)) + ") + " + getBound_Start((*(*lit)->LIS->LB)[i].first);
                } else if ((*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SGE || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_UGE) {
					errs() << "rand() % (" + (getBound_Start((*(*lit)->LIS->LB)[i].first) + " - " + getBound_Start((*(*lit)->LIS->LB)[i].second)) + " + 1) + " + getBound_Start((*(*lit)->LIS->LB)[i].second);
				} else if ((*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_SGT || (*(*lit)->LIS->PREDICATE)[i] == llvm::CmpInst::ICMP_UGT) {
					errs() << "rand() % (" + (getBound_Start((*(*lit)->LIS->LB)[i].first) + " - " + getBound_Start((*(*lit)->LIS->LB)[i].second)) + ") + " + getBound_Start((*(*lit)->LIS->LB)[i].second);
				}
				else {
                    errs() << "\n Error in generating random sample \n";
                }

                errs() << ";\n";
                errs() << space << "    if (" << indvName[(*(*lit)->LIS->IDV)[i]] << "_Start % " <<  getLoopInc((*(*lit)->LIS->INC)[i]) << " != 0) goto SAMPLE; \n";
                if (find(outloops.begin(), outloops.end(), *lit) != outloops.end()) {
                    errs() << space << "    if (" << indvName[(*(*lit)->LIS->IDV)[i]] << "_Start + THREAD_NUM * CHUNK_SIZE > " << getBound_Start((*(*lit)->LIS->LB)[i].second) << ") { goto SAMPLE; }\n";
                }
            }
        }
        
        string idx_string_tmp = "";
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                idx_string_tmp += "std::to_string(";
                idx_string_tmp += indvName[(*(*lit)->LIS->IDV)[i]] + "_Start";
                idx_string_tmp += ") + \"_\" + ";
            }
        }
        idx_string_tmp.pop_back();
        idx_string_tmp.pop_back();
        
        space += "    ";
        errs() << space + "string idx_string = " + idx_string_tmp + ";\n";

/*
        errs() << space + "while ( record.find(idx_string) != record.end() ) {\n";
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                errs() << space + "    " + indvName[(*(*lit)->LIS->IDV)[i]] + "_Start" + " = ";
                errs() << "rand() % (" + (getBound_Start((*(*lit)->LIS->LB)[i].second) + " - " + getBound_Start((*(*lit)->LIS->LB)[i].first)) + ") + " + getBound_Start((*(*lit)->LIS->LB)[i].first);
                errs() << ";\n";
            }
        }
        
        idx_string_tmp = "";
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                idx_string_tmp += "std::to_string(";
                idx_string_tmp += indvName[(*(*lit)->LIS->IDV)[i]] + "_Start";
                idx_string_tmp += ") + \"_\" + ";
            }
        }
        idx_string_tmp.pop_back();
        idx_string_tmp.pop_back();
 
        errs() << space + "    idx_string = " + idx_string_tmp + ";\n";
        errs() << space + "}\n";
*/
        
        errs() << space + "if ( record.find(idx_string) != record.end() ) goto SAMPLE;\n";
        
        errs() << space + "record.insert( idx_string );\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space + "cout << \"[" + refName + to_string(useID) + "]Samples: \" << idx_string << endl;\n";
        errs() << "#endif\n";
#endif
  
        errs() << space + "uint64_t cnt = 0;\n";
        errs() << space + "bool cntStart = false;\n";
        errs() << "\n";
#ifdef PARALLEL
        errs() << space + "/* Variable used to compute thread-local iteration space */\n";
        errs() << space + "auto BLIST = new int[THREAD_NUM][2];\n";
        errs() << space + "int t_Start = 0;\n";
#endif
        
        errs() << space + "/* Generating reuse search code */\n";
        errs() << space + "/* Sampled IDVs " << sampleIDVs.size() << "  */\n";
        for(vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = sampleIDVs.begin(); it != sampleIDVs.end(); ++it) {
            errs() << space + "/* Sampled IDV: " + indvName[(*(*it)->LIS->IDV)[0]] + "  */\n";
        }
        errs() << space + "/* Vector that contains the interleaved iteration, avoid duplicate declaration */\n";
        errs() << space + "vector<vector<int>> nv(THREAD_NUM);\n";
        errs() << space + "int chunk_size, chunk_num, c_Start, ci_Start;\n";
        refRTSearchGen(LoopRefTree, false, true, refName, useID, loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(), sampleIDVs, "    ");
        errs() << "EndSample:\n";
        errs() << space + "s++;\n";
        errs() << space + "}\n";
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        /*
        for (std::map<string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < it->second; i++) {
                errs() << "void ref_" + it->first + std::to_string(i) + "() {\n";
                
//                errs() << "std::cout << \" Start " + it->first + std::to_string(i) + " \";\n " ;
                
                refRTBodyGen(LoopRefTree, it->first, i);
                
//                errs() << "std::cout << \" End " + it->first + std::to_string(i) + " \";\n " ;
                
                errs() << "}\n";
            }
        }
         */
        for (std::map<Instruction*, int>::iterator it = refNumber.begin(), eit = refNumber.end(); it != eit; ++it) {
            errs() << "void ref_" + arrayName[it->first] + std::to_string(it->second) + "() {\n";
            refRTBodyGen(LoopRefTree, arrayName[it->first], it->second);
            errs() << "}\n";
        }
        
        return;
    }
    
    void AccLevelUISamplingCodeGen_ref::mainGen() {
    
        errs() << "int main() {\n";
        
        string space = "    ";
        /*
        for (std::map<string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < (*it).second; i++) {
         */
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "    PAPI_timer_init();\n";
        errs() << "    PAPI_timer_start();\n";
        errs() << "#endif\n";
        errs() << space + "/* " << refNumber.size() << " */\n";
        for (std::map<Instruction*, int>::iterator it = refNumber.begin(), eit = refNumber.end(); it != eit; ++it) {

#ifdef PARALLEL_CXX_THREAD
            /*
                errs() << space + "std::thread t_"+ (*it).first + "_" +std::to_string(i)+ "(";
                errs() << "ref_" + (*it).first + std::to_string(i) + ");\n";
             */
            errs() << space + "std::thread t_"+ arrayName[it->first] + "_" + std::to_string(it->second) + "(";
            errs() << "ref_" + arrayName[it->first] + std::to_string(it->second) + ");\n";
#else
            /*
                errs() << space + "ref_" + (*it).first + std::to_string(i) + "();\n";
             */
            errs() << space + "ref_" + arrayName[it->first] + std::to_string(it->second) + "();\n";
#endif
        /*
            }
         */
        }
        
#ifdef PARALLEL_CXX_THREAD
        /*
        for (std::map<std::string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < (*it).second; i++) {
                    errs() << space + "t_" + (*it).first + "_" + std::to_string(i) + ".join();\n";
            }
        }*/
        for (std::map<Instruction*, int>::iterator it = refNumber.begin(), eit = refNumber.end(); it != eit; ++it) {
            errs() << space + "t_" + arrayName[it->first] + "_" + std::to_string(it->second) + ".join();\n";
        }
#endif

#ifdef UNIFORM_SMOOTHING
        errs() << "    group_uniform_smoothing();\n";
#elif defined(GAUSSIAN_SMOOTHING)
        errs() << "    group_gaussian_smoothing(8.0);\n";   
#endif
        

#ifdef DumpRTMR
#ifdef REFERENCE_GROUP  
        errs() << "    rtMerge();\n";
#endif
        errs() << "    RTtoMR_AET();\n";
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "    PAPI_timer_end();\n";
        errs() << "    PAPI_timer_print();\n";
        errs() << "#endif\n";
        errs() << "    rtDump();\n";
#ifdef REFERENCE_GROUP  
        errs() << "    refRTDump();\n";
#endif
        errs() << "    dumpMR();\n";
#elif defined(DumpRefLease)
        errs() << "    RL_main(0);\n";
#endif
        
        errs() << "    return 0;\n";
        errs() << "}\n";
        
        return;
    }
    
    
    bool AccLevelUISamplingCodeGen_ref::runOnFunction(Function &F) {
    
        errs() << " // Start to generating Static Sampling Code (reference based)\n";
        
        /* reading info from previous passes */
        arrayName = getAnalysis<idxAnalysis::IndexAnalysis>().arrayName;
        arrayExpression = getAnalysis<idxAnalysis::IndexAnalysis>().arrayExpression;
        arrayAccessVariable= getAnalysis<idxAnalysis::IndexAnalysis>().arrayAccessVariable;
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().PTLoopRefTree;
        sampleNum = getAnalysis<sampleNumAnalysis::SampleNumberAnalysis>().sampleNum;

        /* init */
        initArrayName();
        numberRefToSameArray(LoopRefTree);
        numberLoops(LoopRefTree);
        initIndvName(LoopRefTree);
#ifdef PARALLEL
        outloops = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().outMostLoops;
#endif
        /* generate headers */
        headerGen();

        /* generate rtHistoCal function */
        rtHistoGen();

        /* generate subBlkRT function */
        subBlkRTGen();

#if defined(REFERENCE_GROUP)
        rtMergeGen();
#endif

#if defined(UNIFORM_SMOOTHING)
        UniformDistrGen();
#elif defined(GAUSSIAN_SMOOTHING)
        GaussianDistrGen();
#endif
        
#if defined(REFERENCE_GROUP) 
#if defined(GAUSSIAN_SMOOTHING)
        GroupGaussianDistrGen("    ");
#elif defined(UNIFORM_SMOOTHING)
        GroupUniformDistrGen("    ");
#endif
#endif


#ifdef DumpRTMR
        /* generate rtToMR function */
        rtToMRGen();
        
        /* generate rtDump function */
        rtDumpGen();
        
        /* generate mrDump function */
        mrDumpGen();
#elif defined(DumpRefLease)
        /* genearte accessRatioCal function */
        accessRatioCalGen();
        
        /* generate initHitsCosts function */
        initHitsCostsGen();
        
        /* generate getPPUC function */
        getPPUCGen();
        
        /* generate getMaxPPUC function */
        getMaxPPUCGen();
        
        /* generate DumpRI function */
        DumpRIGen();
        
        /* generate RLGen function */
        RLGen();
#endif
        
        /* generate addr cal function*/
        addrCalFuncGenTop(LoopRefTree);

#ifdef SEARCH_REUSE_DIFFERENT_LOOPS
        /* generate  */
        searchReuseDifferentLoopsUpdateFuncGen();
        searchReuseDifferentLoopsCalFuncGen();
#endif
        
        /* generate rtGen */
        refRTGen(LoopRefTree);
        
        /* generate main function */
        mainGen();

        return false;
    }
    
    void AccLevelUISamplingCodeGen_ref::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
        AU.addRequired<loopTreeTransform::ParallelLoopTreeTransform>();
        AU.addRequired<sampleNumAnalysis::SampleNumberAnalysis>();
        return;
    }
    
}
