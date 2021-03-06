#include "ssCodeGen_ref.hpp"


using namespace std;

namespace ssCodeGen_ref {

    char StaticSamplingCodeGen_ref::ID = 0;
    static RegisterPass<StaticSamplingCodeGen_ref> X("ssCodeGen_ref", "static sampling code generating pass (reference based)", false, false);

    StaticSamplingCodeGen_ref::StaticSamplingCodeGen_ref() : FunctionPass(ID) {}

    void StaticSamplingCodeGen_ref::numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree->AA != NULL) {
            /*
            if (refToSameArrayCnt.find(arrayName[LoopRefTree->AA]) == refToSameArrayCnt.end()) {
                refNumber[LoopRefTree->AA] = 0;
                refToSameArrayCnt[arrayName[LoopRefTree->AA]] = 1;
            } else {
                refNumber[LoopRefTree->AA] = refToSameArrayCnt[arrayName[LoopRefTree->AA]];
                refToSameArrayCnt[arrayName[LoopRefTree->AA]] += 1；
            }
            */
            refNumber[LoopRefTree->AA] = refGlobalNumber;
            refGlobalNumber++;
        }
        
        if (LoopRefTree->next != NULL) {
            for ( vector<TreeNodeBase*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(*it);
                if (node != nullptr) {
                    numberRefToSameArray(node);
                }
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree->L != NULL) {
            if (loopNumber.find(LoopRefTree->L) == loopNumber.end()) {
                loopNumber[LoopRefTree->L] = loopNumber.size();
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for ( vector<TreeNodeBase*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(*it);
                if (node != nullptr) {
                    numberLoops(node);
                }
            }
        }
        return;
    }
    
    void StaticSamplingCodeGen_ref::initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree == NULL) {
            return;
        }
        
        if (LoopRefTree->L != NULL) {
            for ( vector<Value *>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvName[*it] = (*it)->getName();
                if ( find(indvName[*it].begin(), indvName[*it].end(), '.') != indvName[*it].end()) {
                    indvName[*it].erase( find(indvName[*it].begin(), indvName[*it].end(), '.'));
                }
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for ( vector<TreeNodeBase*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(*it);
                if (node != nullptr) {
                    initIndvName(node);
                }
                
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::initArrayName() {
        
        for ( map<Instruction*,  string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            it->second.replace( find(it->second.begin(), it->second.end(), '.'),  find(it->second.begin(), it->second.end(), '.') +1, 1, '_');
        }
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree,  vector<string> indvs) {
        
        if (LoopRefTree->L != NULL) {
            for ( vector<Value*>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvs.push_back(indvName[(*it)]);
            }
        }
        
        if (LoopRefTree->AA != NULL) {
			errs() << "/* " + arrayName[LoopRefTree->AA] + " " + arrayExpression[LoopRefTree->AA] + " " +  to_string(refNumber[LoopRefTree->AA])  + " */\n";
            errs() << "int calAddr" + arrayName[LoopRefTree->AA] +  to_string(refNumber[LoopRefTree->AA]) + "( ";
            string arguments = "";
            for ( vector<string>::iterator it = indvs.begin(), eit = indvs.end(); it != eit; ++it) {
                arguments += "int " + (*it) + ", ";
            }
            if (arguments.size() != 0) {
                arguments.pop_back();
                arguments.pop_back();
            }
            errs() << arguments;
            errs() << ") {\n";
#if defined (CLS) && defined (DS)
            errs() << "    int result = (" + arrayExpression[LoopRefTree->AA] + ") * " +  to_string(DS) + " / " +  to_string(CLS) + ";\n";
#else
            errs() << "    int result = " + arrayExpression[LoopRefTree->AA] + ";\n";
#endif
            errs() << "    return result;\n";
            errs() << "}\n";
            return;
        }
        
        if (LoopRefTree->next != NULL) {
            for ( vector<TreeNodeBase*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(*it);
                if (node != nullptr) {
                    addrCalFuncGen(node, indvs);
                }
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        vector<string> indvs;
        addrCalFuncGen(LoopRefTree, indvs);
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::headerGen() {
        
        errs() << "#include <map>\n";
        errs() << "#include <set>\n";
        errs() << "#include <cstdlib>\n";
        errs() << "#include <iostream>\n";
        errs() << "#include <cmath>\n";

#ifdef PARALLEL_CXX_THREAD
        errs() << "#include <thread>\n";
        errs() << "#include <mutex>\n";
#endif
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "#  include <chrono>\n";
        errs() << "#endif\n";
        errs() << "using namespace std;\n";
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "using namespace  chrono;\n";
        errs() << "#endif\n";
        
#ifdef PARALLEL_CXX_THREAD
        errs() << " mutex mtx;\n";
#endif
        
        return;
    }

#ifdef DumpRTMR
    void StaticSamplingCodeGen_ref::rtHistoGen() {
        
        errs() << " map<uint64_t, double> RT;\n";
        errs() << " map<uint64_t, double> MR;\n";
        errs() << "void rtHistoCal( int rt) {\n";
    #ifdef PARALLEL_OMP
        errs() << "    #pragma omp critical\n";
        errs() << "    {";
        errs() << "        if (RT.find(rt) == RT.end()) { \n";
        errs() << "            RT[rt] = 1;\n";
        errs() << "        } else {\n";
        errs() << "            RT[rt] += 1;\n";
        errs() << "        }\n";
        errs() << "    }";
    #elif defined(PARALLEL_CXX_THREAD)
        errs() << "     unique_lock< mutex> lck (mtx, defer_lock);\n";
        errs() << "    lck.lock();\n";
        errs() << "    if (RT.find(rt) == RT.end()) { \n";
        errs() << "        RT[rt] = 1;\n";
        errs() << "    } else {\n";
        errs() << "        RT[rt] += 1;\n";
        errs() << "    }\n";
        errs() << "    lck.unlock();\n";
    #else
        errs() << "    if (RT.find(rt) == RT.end()) { \n";
        errs() << "        RT[rt] = 1;\n";
        errs() << "    } else {\n";
        errs() << "        RT[rt] += 1;\n";
        errs() << "    }\n";
    #endif
        errs() << "    return;\n";
        errs() << "}\n";
    }
#elif defined(DumpRefLease)
    void StaticSamplingCodeGen_ref::rtHistoGen() {
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
    
    
#ifdef DumpRTMR
    void StaticSamplingCodeGen_ref::rtDumpGen() {
        
        errs() << "void rtDump() {\n";
        errs() << "    cout << \"Start to dump reuse time histogram\\n\";\n";
        errs() << "    for (map<uint64_t, double>::iterator it = RT.begin(), eit = RT.end(); it != eit; ++it) {\n";
        errs() << "        cout << it->first << \" \" << it->second << \"\\n\";\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::rtToMRGen() {
        
        errs() << "void RTtoMR_AET() {\n";
        
        string space = "    ";
        errs() << space + " map<uint64_t, double> P;\n";
        errs() << space + "double total_num_RT = 0;\n";
        errs() << space + "uint64_t max_RT = 0;\n";
        
        errs() << space + "for ( map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {\n";
        errs() << space + "    total_num_RT += it->second;\n";
        errs() << space + "    if (max_RT < it->first) {\n";
        errs() << space + "        max_RT = it->first;\n";
        errs() << space + "    }\n";
        errs() << space + "}\n";
        
        errs() << space + "double accumulate_num_RT = 0;\n";
        errs() << space + "for ( map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {\n";
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
    
    void StaticSamplingCodeGen_ref::mrDumpGen() {
        
        errs() << "void dumpMR() {\n";
        
        errs() << "    cout << \"miss ratio\" << endl;\n";
        
        errs() << "     map<uint64_t, double>::iterator it1 = MR.begin();\n";
        errs() << "     map<uint64_t, double>::iterator it2 = MR.begin();\n";
        
        errs() << "    while(it1 != MR.end()) {\n";
        errs() << "        while(1) {\n";
        errs() << "             map<uint64_t, double>::iterator it3 = it2;\n";
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
        errs() << "        cout << it1->first << \" \" << it1->second << endl;\n";
        errs() << "        if (it1 != it2) {\n";
        errs() << "            cout << it2->first << \" \" << it2->second << endl;\n";
        errs() << "        }\n";
        errs() << "        it1 = ++it2;\n";
        errs() << "        it2 = it1;\n";
        errs() << "    }\n";
        
        errs() << "    return;\n";
        errs() << "}\n";
        
        return;
    }
#elif defined(DumpRefLease)
    void StaticSamplingCodeGen_ref::accessRatioCalGen() {
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
    
    void StaticSamplingCodeGen_ref::initHitsCostsGen() {
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
    void StaticSamplingCodeGen_ref::getPPUCGen() {
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
    void StaticSamplingCodeGen_ref::getMaxPPUCGen() {
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
	void StaticSamplingCodeGen_ref::DumpRIGen() {
        errs() << "void dumpRI() {\n";
        errs() << "    uint64_t total_number_of_ri = 0;\n";
        errs() << "    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {\n";
        errs() << "         set<uint64_t> riset;\n";
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
    void StaticSamplingCodeGen_ref::RLGen() {
        errs() << "void RL_main(uint64_t CacheSize) {\n";
        errs() << "    initHitsCosts();\n";
        errs() << "    accessRatioCal();\n";
        errs() << "    double totalCost = 0;\n";
        errs() << "    double totalHitRatio = 0;\n";
        errs() << "    double targetCost = CacheSize;\n";
        errs() << "#ifdef DEBUG\n";
        errs() << "    dumpRI();\n";
        errs() << "#endif\n";
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

    string StaticSamplingCodeGen_ref::getBound(Value *bound) {
        
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
    
    string StaticSamplingCodeGen_ref::getBound_Start(Value *bound) {
        
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
    
    vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> StaticSamplingCodeGen_ref::findLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, string refName, int useID, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops) {
        
        if (LoopRefTree->L != NULL) {
            loops.push_back(LoopRefTree);
        }
        if (LoopRefTree->AA != NULL) {
            if (arrayName[LoopRefTree->AA] == refName && refNumber[LoopRefTree->AA] == useID) {
                return loops;
            } else {
                return  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>();
            }
        }
        
         vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopRes;
        
        if (LoopRefTree->next != NULL) {
            for ( vector<TreeNodeBase*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(*it);
                if (node != nullptr) {
                    vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopTmp = findLoops(node, refName, useID, loops);
                    if (loopTmp.size() != 0) {
                        loopRes = loopTmp;
                    }
                }
            }
        }
        return loopRes;
    }
    
    /* Search result reuse (Different loop) */
    /* a is coefficient vector */
    /* x is the input matrix */
    /* y is the input result */
    /* length */
    void StaticSamplingCodeGen_ref::searchReuseDifferentLoopsUpdateFuncGen() {
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
    
    void StaticSamplingCodeGen_ref::searchReuseDifferentLoopsCalFuncGen() {
        
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
    
    bool StaticSamplingCodeGen_ref::searchReuseDifferentLoopsInitGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, bool GenFlag,  string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> currentLoops, string space) {
        
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
                    errs() << space + "int* prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = new int[" +  to_string(numRecord) +"];\n";
                    errs() << space + " fill(prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + ", prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " + " +  to_string(numRecord) + ", -1);\n";
                    
                    errs() << space + "int** x_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = new int*["+  to_string(numRecord) +"];\n";
                    errs() << space + "float* a_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = new float["+  to_string(numRecord) +"];\n";
                    errs() << space + "int* curr_v_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = new int["+  to_string(numRecord) +"];\n";
                    int cntTmp = 0;
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        errs() << space + "int* prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = new int[" +  to_string(numRecord) +"];\n";
                        errs() << space + " fill(prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + ", prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " + " +   to_string(numRecord) + ", -1);\n";
                        errs() << space + "x_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[" +  to_string(cntTmp) + "] = prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + ";\n";
                        cntTmp++;
                    }
                    
                    errs() << space + "int * prev_one_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = new int[" +  to_string(numRecord) +"];\n";
                    errs() << space + " fill(prev_one_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + ", prev_one_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " + " +  to_string(numRecord) + ", 1);\n";
                    errs() << space + "x_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[" +  to_string(cntTmp) + "] = prev_one_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + ";\n";
                    
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        errs() << space + "int prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[" +  to_string(numRecord) +"];\n";
                        errs() << space + "float* a_" + indvName[(*(*it)->LIS->IDV)[0]] + "_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = new float["+  to_string(numRecord) +"];\n";
                        errs() << space + " fill(prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName +  to_string(refNumber[LoopRefTree->AA]) + ", prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " + " +   to_string(numRecord) + ", -1);\n";
                    }
                    errs() << space + "int prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "_idx = 0;\n";
                }
            }
        }
        
        if (LoopRefTree->next != NULL) {
            if (LoopRefTree == loops.back()) {
                return GenFlag;
            }
            for ( vector<TreeNodeBase*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops_New = currentLoops;
                if (LoopRefTree->L != NULL) {
                    currentLoops_New.push_back(LoopRefTree);
                }
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(*it);
                if (node == nullptr) { continue; }
                GenFlag = searchReuseDifferentLoopsInitGen(node, GenFlag, refName, useID, loops, currentLoops_New, space);
            }
        }
        
        return GenFlag;
    }
    
    bool StaticSamplingCodeGen_ref::searchReuseDifferentLoopsBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag,  string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> currentLoops, string space) {
        
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
                     vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator firstIndv = loops.begin();
                    for (unsigned long i = 0; i < currentLoops.size(); i++) {
                        errs() << "prev_" + indvName[(*(*firstIndv)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[" +  to_string(i) + "] != -1";
                        if (i + 1 != currentLoops.size()) {
                            errs() << " && ";
                        }
                    }
                    errs() << ") {\n";

                    /* generate assignment statements to init curr_v_ */
                    for (unsigned long i = 0; i < loops.size(); i++) {
                        errs() << "curr_v_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[" +  to_string(i) + "] = " + indvName[(*(loops[i])->LIS->IDV)[0]] + "_Start;\n";
                    }
                    errs() << "curr_v_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[" +  to_string(loops.size()) + "] = 1;\n";
                    
                    /* generate if statement to check prediction access the same location or not */
                    errs() << space + "    if ( calAddr" + refName +  to_string(refNumber[LoopRefTree->AA]) + "( ";
                    for (unsigned long i = 0; i < currentLoops.size(); i++) {
                        errs() << "calWithCoefficient( ";
                        for (unsigned long j = 0; j < loops.size(); j++) {
                            
                        }
                        errs() << "a_" + indvName[(*(currentLoops[i])->LIS->IDV)[0]] + "_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                        errs() << " , ";
                        errs() << "curr_v_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                        errs() << " , ";
                        errs() <<  to_string(loops.size());
                        errs() << ")";
                        if (i + 1 != currentLoops.size()) {
                            errs() << ",";
                        }
                    }
                    errs() << + ") == ";

                    errs() << "calAddr" + refName +  to_string(useID) + "( ";
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
                    errs() << "a_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                    errs() << " , ";
                    errs() << "curr_v_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                    errs() << " , ";
                    errs() <<  to_string(loops.size());
                    errs() << ");\n";
#endif
                    
                    /* predict cnt and accumulate */
                    errs() << space + "        rtHistoCal(";
                    errs() << " calWithCoefficient(";
                    errs() << "a_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                    errs() << " , ";
                    errs() << "curr_v_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                    errs() << " , ";
                    errs() <<  to_string(loops.size());
                    errs() << ")";
#ifdef DumpRefLease
                    errs() << ", " +  to_string(refNumber[LoopRefTree->AA]);
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
            for ( vector<TreeNodeBase*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops_New = currentLoops;
                if (LoopRefTree->L != NULL) {
                    currentLoops_New.push_back(LoopRefTree);
                }
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(*it);
                if (node != nullptr) {
                    GenFlag = searchReuseDifferentLoopsBodyGen(node, GenFlag, refName, useID, loops, currentLoops_New, space);
                }
            }
        }
        
        return GenFlag;
    }
    

    /* Search result reuse (Same loop) */
    void StaticSamplingCodeGen_ref::searchReuseSameLoopInitGen( string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, string space) {
        
        //errs() << "Search result reuse: (ref name " << refName << " ) ( ID " << useID << " ) ( numOfLoops " << loops.size() << " )\n";

        for ( vector<TreeNodeBase*>::iterator ait = loops.back()->next->begin(), eit = loops.back()->next->end(); ait != eit; ++ait) {
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(*ait);
            if (node == nullptr) { continue; }
            if (node->AA != NULL) {
                if (arrayName[node->AA] == refName) {
                    errs() << space + "uint64_t prev_cnt_" + refName +  to_string(refNumber[node->AA]) + " = -1;\n";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        errs() << space + "uint64_t prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[node->AA]) + " = -1;\n";
                        errs() << space + "uint64_t prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName +  to_string(refNumber[node->AA]) + " = -1;\n";
                    }
                }
            }
        }
    
        return;
    }
    
    void StaticSamplingCodeGen_ref::searchReuseSameLoopBodyGen( string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops,  string space) {
        
        for ( vector<TreeNodeBase*>::iterator ait = loops.back()->next->begin(), eit = loops.back()->next->end(); ait != eit; ++ait) {
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(*ait);
            if (node == nullptr) { continue; }
            if (node->AA != NULL) {
                if (arrayName[node->AA] == refName) {
                    errs() << space + "if ( prev_cnt_" + refName +  to_string(refNumber[node->AA]) + " != -1) {\n";
                    errs() << space + "    if ( calAddr" + refName +  to_string(refNumber[node->AA]) + "( ";
                    string tmp = "";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            
                            tmp += indvName[(*(*it)->LIS->IDV)[i]] + "_Start";
                            tmp += " - prev_" + indvName[(*(*it)->LIS->IDV)[i]] + "_Start_" + refName +  to_string(refNumber[node->AA]);
                            tmp += " + prev_" + indvName[(*(*it)->LIS->IDV)[i]] + "_End_" + refName +  to_string(refNumber[node->AA]);
                            tmp += ", ";
                        }
                    }
                    if (loops.size() != 0) {
                        tmp.pop_back();
                        tmp.pop_back();
                    }
                    
                    errs() << tmp + ")";
                    errs() << " == ";
                    
                    errs() << "calAddr" + refName +  to_string(useID);
                    errs() << "(";
                    
                    tmp = "";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
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
                    errs() << space + "        pred = " + "prev_cnt_" + refName +  to_string(refNumber[node->AA]) + ";\n";
#endif
                    
                    errs() << space + "        rtHistoCal(prev_cnt_" + refName +  to_string(refNumber[node->AA]);
#ifdef DumpRefLease
                    errs() << ", " +  to_string(refNumber[node->AA]);
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
    bool StaticSamplingCodeGen_ref::refRTSearchGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag, string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, string space) {
        
        if (loops.size() != 0 && GenFlag == false) {
            if (LoopRefTree == loops[0]) {
                GenFlag = true;
            }
        }
        
        if (GenFlag == true) {
            
            /* Generate loop */
            if (LoopRefTree->L != NULL) {
                
                errs() << space + "{\n";
                
                string loopNum =  to_string(loopNumber[LoopRefTree->L]);
                if ( find(loops.begin(), loops.end(), LoopRefTree) != loops.end()) {
                    if (LoopRefTree == loops[0]) {
                        errs() << space + "int " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = ";
                        errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] + "_Start;\n";
                    } else {
                        errs() << space + "int " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = ";
                        errs() << getBound((*LoopRefTree->LIS->LB)[0].first) + ";\n";
                        errs() << space + "if ( ";
                        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                            if (LoopRefTree == *it) {
                                break;
                            }
                            if (it != loops.begin()) {
                                errs() << " && ";
                            }
                            errs() << indvName[(*(*it)->LIS->IDV)[0]] << " == " + indvName[(*(*it)->LIS->IDV)[0]] + "_Start";
                        }
                        errs() << " ) {\n";
                        errs() << space + "    " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = ";
                        errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] + "_Start;\n";
                        errs() << space + "}\n";
                    }
                } else {
                    errs() << space + "int " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = ";
                    errs() << getBound((*LoopRefTree->LIS->LB)[0].first) + ";\n";
                }
                
                errs() << space + "for ( int " + indvName[(*LoopRefTree->LIS->IDV)[0]];
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
                
                errs() << getBound((*LoopRefTree->LIS->LB)[0].second);
                errs() << "; ";
                
                /* need to take stride into consideration */
                if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULT) {
                    
                    errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] + "++";
                    
                } else if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGT) {
                    
                    errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] + "--";
                    
                } else {
                    errs() << "\n Error in geting stride \n";
                }
                
                
                errs() << ") {\n";
            }
            
            /* Generate addr checking instruction */
            if (LoopRefTree->AA != NULL) {
                
                if (arrayName[LoopRefTree->AA] == refName) {
                
                    
                    errs() << space + "if (cntStart == true) {\n";
                    errs() << space + "    cnt++;\n";
                    
                    errs() << space + "    if ( calAddr" + refName +  to_string(refNumber[LoopRefTree->AA]) + "( ";
                    string tmp = "";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            tmp += indvName[(*(*it)->LIS->IDV)[i]];
                            tmp += ", ";
                        }
                    }
                    if (currentLoops.size() != 0) {
                        tmp.pop_back();
                        tmp.pop_back();
                    }
                    
                    errs() << tmp + ")";
                    errs() << " == ";
                    
                    errs() << "calAddr" + refName +  to_string(useID);
                    errs() << "(";
                    
                    tmp = "";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
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
                    
                    errs() << space + "        rtHistoCal(cnt";
#ifdef DumpRefLease
                    errs() << ", " +  to_string(useID);
#endif
                    errs() << ");\n";
                    
#ifdef PROFILE_SEARCH_REUSE
                    errs() << space + "        actural = cnt;\n";
#endif

#ifdef SEARCH_REUSE_SAME_LOOP
                    if ( find(loops.back()->next->begin(), loops.back()->next->end(), LoopRefTree) != loops.back()->next->end()) {
#ifdef PROFILE_SEARCH_REUSE
                        errs() << space + "        reuse_in_same_loop++;\n";
#endif
                        errs() << space + "        prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = cnt;\n";
                        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                            errs() << space + "        prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = " + indvName[(*(*it)->LIS->IDV)[0]] + "_Start" + ";\n";
                            errs() << space + "        prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName +  to_string(refNumber[LoopRefTree->AA]) + " = " + indvName[(*(*it)->LIS->IDV)[0]] + ";\n";
                        }
                    }
#endif

                    
#ifdef SEARCH_REUSE_DIFFERENT_LOOPS
//                    if ( find(loops.back()->next->begin(), loops.back()->next->end(), LoopRefTree) == loops.back()->next->end())  {
                    if (loops.front() != currentLoops.front()) {
#ifdef PROFILE_SEARCH_REUSE
                        errs() << space + "        reuse_in_diff_loops++;\n";
#endif
                        errs() << space + "        prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "_idx" + "] = cnt;\n";
                        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                            errs() << space + "        prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "_idx" + "] = " + indvName[(*(*it)->LIS->IDV)[0]] + "_Start" + ";\n";
                        }
                        errs() << space + "prev_one_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "_idx" + "] = 1;\n";
                        
                        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                            errs() << space + "        prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "_idx" + "] = " + indvName[(*(*it)->LIS->IDV)[0]] + ";\n";
                        }
                        
                        errs() << space + "        prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "_idx = (prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "_idx + 1) % " +  to_string(loops.size() + 1) + " ;\n";
                        
                        errs() << space + "        if ( ";
                         vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator firstIndv = loops.begin();
                        for (unsigned long i = 0; i < currentLoops.size(); i++) {
                            errs() << "prev_" + indvName[(*(*firstIndv)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[LoopRefTree->AA]) + "[" +  to_string(i) + "] != -1";
                            if (i + 1 != currentLoops.size()) {
                                errs() << " && ";
                            }
                        }
                        
                        errs() << ") {\n";
                        for (unsigned long i = 0; i < currentLoops.size(); i++) {
                            errs() << space + "           updateCoefficient(";
                            errs() << "a_" + indvName[(*(currentLoops[i])->LIS->IDV)[0]] + "_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                            errs() << ", x_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                            errs() << ", " +  to_string(loops.size() + 1);
                            errs() << ", prev_" + indvName[(*(currentLoops[i])->LIS->IDV)[0]] + "_End_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                            errs() << ");\n";
                        }
                        errs() << space + "           updateCoefficient(";
                        errs() << "a_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                        errs() << ", x_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                        errs() << ", " +  to_string(loops.size() + 1);
                        errs() << ", prev_cnt_" + refName +  to_string(refNumber[LoopRefTree->AA]);
                        errs() << ");\n";
                        errs() << space + "        }\n";
                    }
#endif
                    errs() << space + "        goto EndSample;\n";
                    
                    errs() << space + "    }\n";
                    errs() << space + "}\n";
                    
                    if (useID == refNumber[LoopRefTree->AA]) {
                        errs() << space + "cntStart = true;\n";
                    }
                } else {
                    errs() << space + "if (cntStart == true) cnt++;\n";
                }
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for ( vector<TreeNodeBase*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops_New = currentLoops;
                if (LoopRefTree->L != NULL) {
                    currentLoops_New.push_back(LoopRefTree);
                }
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(*it);
                if (node != nullptr) {
                    GenFlag = refRTSearchGen(node, GenFlag, refName, useID, loops, currentLoops_New, space + "    ");
                }
            }
            if (LoopRefTree->L != NULL && GenFlag == true) {
                errs() << space + "}\n";
                errs() << space + "}\n";
            }
        }
        
        return GenFlag;
    }


    void StaticSamplingCodeGen_ref::refRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, string refName, int useID) {
        
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops = findLoops(LoopRefTree, refName, useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>());
        
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
        
#ifdef SEARCH_REUSE_SAME_LOOP
        errs() << "    /* Generating search reuse init code (same loop) */\n";
        if (loops.size() != 0) {
            searchReuseSameLoopInitGen(refName, useID, loops, "    ");
        }
#endif

#ifdef SEARCH_REUSE_DIFFERENT_LOOPS
        errs() << "    /* Generating search reuse init code (different loops) */\n";
        if (loops.size() != 0) {
            searchReuseDifferentLoopsInitGen(LoopRefTree, false, refName, useID, loops,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(), "    ");
        }
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
            	errs() <<  to_string(sampleNum[loops.back()]);
        	} else {
            	errs() << "ERROR in finding bounds\n";
        	}
		}
        errs() << ";) {\n";
        
        if (loops.size() == 0) {
            errs() << space + "/* Generating reuse search code */\n";
            errs() << "\n";
            refRTSearchGen(LoopRefTree, false, refName, useID, loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(), "    ");
            errs() << space + "s++;\n";
            errs() << space + "}\n";
            return;
        }

        errs() << "SAMPLE:\n";
        
        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
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

                
            }
        }
        
        string idx_string_tmp = "";
        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                idx_string_tmp += " to_string(";
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
        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                errs() << space + "    " + indvName[(*(*lit)->LIS->IDV)[i]] + "_Start" + " = ";
                errs() << "rand() % (" + (getBound_Start((*(*lit)->LIS->LB)[i].second) + " - " + getBound_Start((*(*lit)->LIS->LB)[i].first)) + ") + " + getBound_Start((*(*lit)->LIS->LB)[i].first);
                errs() << ";\n";
            }
        }
        
        idx_string_tmp = "";
        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                idx_string_tmp += " to_string(";
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
#endif
  
        errs() << space + "uint64_t cnt = 0;\n";
        errs() << space + "bool cntStart = false;\n";
        errs() << "\n";
        
        
#ifdef SEARCH_REUSE_SAME_LOOP
        errs() << "        /* Generating search reuse body code (use reuse are in the same loop) */\n";
        if (loops.size() != 0) {
            searchReuseSameLoopBodyGen(refName, useID, loops, "        ");
        }
#endif
        
#ifdef SEARCH_REUSE_DIFFERENT_LOOPS
        errs() << "        /* Generating search reuse body code (use reuse are in different loop) */\n";
        if (loops.size() != 0) {
            searchReuseDifferentLoopsBodyGen(LoopRefTree, false, refName, useID, loops,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(), "        ");
        }
#endif
        
        
#ifdef PROFILE_SEARCH_REUSE
        errs() << "        /* Finished search reuse */\n";
        errs() << "PROFILE:\n";
        errs() << "sample_num++;\n";
#endif
        
        errs() << space + "/* Generating reuse search code */\n";
        errs() << "\n";
        refRTSearchGen(LoopRefTree, false, refName, useID, loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(), "    ");
        
#ifdef PROFILE_SEARCH_REUSE
        errs() << space + "no_reuse_num++;\n";
        errs() << space + "if (actural != 0 && actural != pred) mis_pred_num++;\n";
        errs() << space + "actural = 0;\n";
#endif
        
        errs() << "EndSample:\n";
        errs() << space + "s++;\n";
        errs() << space + "}\n";

#ifdef PROFILE_SEARCH_REUSE
        // SN: sample number, NSR: Num of samples have reuse, RSL: reuses from same loop, RDL: reuses from diff loops, PN: Predict num, MPN: Miss Predict num
        // PSL: predicted in same loop, PDL: predicted in diff loops
        errs() << " cout << \"SN \" << sample_num << \" \";\n";
        errs() << " cout << \"NSR \" << sample_num - no_reuse_num << \" \";\n";
        errs() << " cout << \"RSL \" << reuse_in_same_loop << \" \";\n";
        errs() << " cout << \"RDL \" << reuse_in_diff_loops << \" \";\n";
        errs() << " cout << \"PN \" << pred_num << \" \";\n";
        errs() << " cout << \"PSL \" << pred_sl_num << \" \";\n";
        errs() << " cout << \"PDL \" << pred_dl_num << \" \";\n";
        errs() << " cout << \"MPN \" << mis_pred_num <<  endl;\n";
#endif
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        /*
        for ( map<string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < it->second; i++) {
                errs() << "void ref_" + it->first +  to_string(i) + "() {\n";
                
//                errs() << " cout << \" Start " + it->first +  to_string(i) + " \";\n " ;
                
                refRTBodyGen(LoopRefTree, it->first, i);
                
//                errs() << " cout << \" End " + it->first +  to_string(i) + " \";\n " ;
                
                errs() << "}\n";
            }
        }
         */
        for ( map<Instruction*, int>::iterator it = refNumber.begin(), eit = refNumber.end(); it != eit; ++it) {
            errs() << "void ref_" + arrayName[it->first] +  to_string(it->second) + "() {\n";
            refRTBodyGen(LoopRefTree, arrayName[it->first], it->second);
            errs() << "}\n";
        }
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::mainGen() {
    
        errs() << "int main() {\n";
        
        string space = "    ";
        /*
        for ( map<string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < (*it).second; i++) {
         */
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "    // Get starting timepoint\n";
        errs() << "    auto start = high_resolution_clock::now();\n" ;
        errs() << "#endif\n";
        for ( map<Instruction*, int>::iterator it = refNumber.begin(), eit = refNumber.end(); it != eit; ++it) {

#ifdef PARALLEL_CXX_THREAD
            /*
                errs() << space + " thread t_"+ (*it).first + "_" + to_string(i)+ "(";
                errs() << "ref_" + (*it).first +  to_string(i) + ");\n";
             */
            errs() << space + " thread t_"+ arrayName[it->first] + "_" +  to_string(it->second) + "(";
            errs() << "ref_" + arrayName[it->first] +  to_string(it->second) + ");\n";
#else
            /*
                errs() << space + "ref_" + (*it).first +  to_string(i) + "();\n";
             */
            errs() << space + "ref_" + arrayName[it->first] +  to_string(it->second) + "();\n";
#endif
        /*
            }
         */
        }
        
#ifdef PARALLEL_CXX_THREAD
        /*
        for ( map< string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < (*it).second; i++) {
                    errs() << space + "t_" + (*it).first + "_" +  to_string(i) + ".join();\n";
            }
        }*/
        for ( map<Instruction*, int>::iterator it = refNumber.begin(), eit = refNumber.end(); it != eit; ++it) {
            errs() << space + "t_" + arrayName[it->first] + "_" +  to_string(it->second) + ".join();\n";
        }
#endif
        

#ifdef DumpRTMR
        errs() << "    rtDump();\n";
        errs() << "    RTtoMR_AET();\n";
        errs() << "    dumpMR();\n";
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "    // Get ending timepoint\n"; 
        errs() << "    auto stop = high_resolution_clock::now(); \n";
  
        errs() << "    // Get duration. Substart timepoints to\n";
        errs() << "    // get durarion. To cast it to proper unit\n"; 
        errs() << "    // use duration cast method\n";
        errs() << "    auto duration = duration_cast<microseconds>(stop - start);\n ";
        errs() << "    cout << \"Time taken by SPS:  \" << duration.count() << endl; \n";
        errs() << "#endif\n";
#elif defined(DumpRefLease)
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "// Get ending timepoint\n"; 
        errs() << "    auto stop = high_resolution_clock::now(); \n";
  
        errs() << "    // Get duration. Substart timepoints to\n";
        errs() << "    // get durarion. To cast it to proper unit\n"; 
        errs() << "    // use duration cast method\n";
        errs() << "    auto duration = duration_cast<microseconds>(stop - start);\n ";
        errs() << "    cout << \"Time taken by SPS:  \" << duration.count() << endl; \n";
        errs() << "#endif\n";
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "    // Get starting timepoint\n";
        errs() << "    start = high_resolution_clock::now();\n" ;
        errs() << "#endif\n";
        errs() << "    RL_main(0);\n";
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "// Get ending timepoint\n"; 
        errs() << "    stop = high_resolution_clock::now(); \n";
  
        errs() << "    // Get duration. Substart timepoints to\n";
        errs() << "    // get durarion. To cast it to proper unit\n"; 
        errs() << "    // use duration cast method\n";
        errs() << "    duration = duration_cast<microseconds>(stop - start);\n ";
        errs() << "    cout << \"Time taken by CARL:  \" << duration.count() << endl; \n";
        errs() << "#endif\n";
#elif
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "// Get ending timepoint\n"; 
        errs() << "    auto stop = high_resolution_clock::now(); \n";
  
        errs() << "    // Get duration. Substart timepoints to\n";
        errs() << "    // get durarion. To cast it to proper unit\n"; 
        errs() << "    // use duration cast method\n";
        errs() << "    auto duration = duration_cast<microseconds>(stop - start);\n ";
        errs() << "    cout << \"Time taken by SPS:  \" << duration.count() << endl; \n";
        errs() << "#endif\n";
#endif
        
        errs() << "    return 0;\n";
        errs() << "}\n";
        
        return;
    }
    
    
    bool StaticSamplingCodeGen_ref::runOnFunction(Function &F) {
    
        errs() << " // Start to generating Static Sampling Code (reference based)\n";
        
        /* reading info from previous passes */
        arrayName = getAnalysis<idxAnalysis::IndexAnalysis>().arrayName;
        arrayExpression = getAnalysis<idxAnalysis::IndexAnalysis>().arrayExpression;
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree;
        sampleNum = getAnalysis<sampleNumAnalysis::SampleNumberAnalysis>().sampleNum;

        /* init */
        initArrayName();
        numberRefToSameArray(LoopRefTree);
        numberLoops(LoopRefTree);
        initIndvName(LoopRefTree);

        /* generate headers */
        headerGen();

        /* generate rtHistoCal function */
        rtHistoGen();
        
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
    
    void StaticSamplingCodeGen_ref::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
        AU.addRequired<sampleNumAnalysis::SampleNumberAnalysis>();
        return;
    }
    
}
