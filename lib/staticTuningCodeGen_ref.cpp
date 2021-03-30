#include "staticTuningCodeGen_ref.hpp"


using namespace std;

namespace staticTuningCodeGen_ref {

    char TuningCodeGen_ref::ID = 0;
    static RegisterPass<TuningCodeGen_ref> X("staticTuningCodeGen_ref", "locality-centric thread number tuning", false, false);

    TuningCodeGen_ref::TuningCodeGen_ref() : FunctionPass(ID) {}

    void TuningCodeGen_ref::numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
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
            for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                if (!(*it)->isThreadNode) {
                    numberRefToSameArray(*it);
                }
            }
        }
        
        return;
    }
    
    void TuningCodeGen_ref::numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree->L != NULL) {
            if (loopNumber.find(LoopRefTree->L) == loopNumber.end()) {
                loopNumber[LoopRefTree->L] = loopNumber.size();
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                numberLoops(*it);
            }
        }
        return;
    }
    
    void TuningCodeGen_ref::initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree == NULL) {
            return;
        }
        
        if (LoopRefTree->L != NULL) {
            for ( vector<Value *>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvName[*it] = (*it)->getName().str();
                if ( find(indvName[*it].begin(), indvName[*it].end(), '.') != indvName[*it].end()) {
                    indvName[*it].erase( find(indvName[*it].begin(), indvName[*it].end(), '.'));
                }
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                initIndvName(*it);
            }
        }
        
        return;
    }
    
    void TuningCodeGen_ref::initArrayName() {
        
        for ( map<Instruction*,  string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            it->second.replace( find(it->second.begin(), it->second.end(), '.'),  find(it->second.begin(), it->second.end(), '.') +1, 1, '_');
        }
        
        return;
    }
    
    void TuningCodeGen_ref::addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree,  vector<string> indvs) {
        
        if (LoopRefTree->L != NULL) {
            for ( vector<Value*>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvs.push_back(indvName[(*it)]);
            }
        }
        
        if (LoopRefTree->AA != NULL) {
            defineRefs(LoopRefTree->AA);
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
            for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                addrCalFuncGen(*it, indvs);
            }
        }
        
        return;
    }
    
    void TuningCodeGen_ref::addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        vector<string> indvs;
        addrCalFuncGen(LoopRefTree, indvs);
        
        return;
    }
    
    void TuningCodeGen_ref::headerGen() {
        
        errs() << "#include <cstdlib>\n";
        errs() << "#include <cmath>\n";
        errs() << "#include <functional>\n";
        errs() << "#include <iostream>\n";
        errs() << "#include <map>\n";
#ifdef PARALLEL_CXX_THREAD
        errs() << "#include <mutex>\n";
#endif
        errs() << "#include <set>\n";
        errs() << "#include <thread>\n";
        errs() << "#include <vector>\n";
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "#  include \"papi_timer.h\"\n";
        errs() << "#endif\n";
        errs() << "using namespace std;\n";
        errs() << "using namespace placeholders;\n";
        // errs() << "#ifdef PAPI_TIMER\n";
        // errs() << "using namespace  chrono;\n";
        // errs() << "#endif\n";
        
#ifdef PARALLEL_CXX_THREAD
        errs() << "mutex mtx;\n";
#endif
        errs() << "map<int, double> per_threadcnt_expected_rt;\n";
        errs() << "map<int, map<uint64_t, double>> per_tcount_RT;\n";
        errs() << "map<int, map<uint64_t, double>> privateRT;\n";
        errs() << "map<int, double> per_tcount_MR;\n";
        errs() << "map<uint64_t, double> RT;\n";
        errs() << "map<uint64_t, double> MR;\n";
        errs() << "uint64_t sample_sum = 0;\n";
#ifdef REFERENCE_GROUP
        errs() << "map<string, map<uint64_t, double>> refRT;\n";
#endif
        errs() << "double total_private_stat = 0.0;\n";
        errs() << "double total_reuse = 0.0;\n";
        errs() << "double total_neighbor = 0.0;\n";
        errs() << "double total_scale = 0.0;\n";
        errs() << "double total_fold = 0.0;\n";
        errs() << "double total_interchunk = 0.0;\n";
        errs() << "double share_reuse = 0.0;\n";

#ifdef SMOOTHING
        errs() << "enum class ReuseType {\n";
        errs() << "    INTER_CHUNK,\n";
        errs() << "    NEIGHBOR,\n";
        errs() << "    SCALE,\n";
        errs() << "    FOLD_SRC_SINK,\n";
        errs() << "    FOLD_SINK_SRC\n";
        errs() << "};\n";
#endif
        return;
    }

    void TuningCodeGen_ref::rtHistoGen() {

        errs() << "void privateRTHistoCal(int tid, map<int, map<uint64_t, double>> &rth, uint64_t rt, double val) {\n";
        errs() << "    if (rth.find(tid) == rth.end()) { \n";
        errs() << "        rth[tid][rt] = val;\n";
        errs() << "    } else {\n";
        errs() << "        rth[tid][rt] += val;\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";

        errs() << "void rtHistoCal( map<uint64_t, double> &rth, uint64_t rt, double val ) {\n";
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
        errs() << "     unique_lock< mutex> lck (mtx, defer_lock);\n";
        errs() << "    lck.lock();\n";
        errs() << "    if (rth.find(rt) == rth.end()) { \n";
        errs() << "        rth[rt] = val;\n";
        errs() << "    } else {\n";
        errs() << "        rth[rt] += val;\n";
        errs() << "    }\n";
        errs() << "    lck.unlock();\n";
#else
        errs() << "    if (rth.find(rt) == rth.end()) { \n";
        errs() << "        rth[rt] = val;\n";
        errs() << "    } else {\n";
        errs() << "        rth[rt] += val;\n";
        errs() << "    }\n";
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

     /* Generate the function to calculate the bins */
    void TuningCodeGen_ref::subBlkRTGen() {
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

    void TuningCodeGen_ref::mrDumpGen() {
        
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
    
    void TuningCodeGen_ref::rtToMRGen() {
        string space = "    ";
        errs() << "void RTtoMR_AET() {\n";
        errs() << space + "map<uint64_t, double> P;\n";
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

        errs() << "double RTtoMR_AET_C(map<uint64_t, double> &rt, uint64_t cache_size) {\n";
        errs() << space + "map<uint64_t, double> P;\n";
        errs() << space + "double total_num_RT = 0;\n";
        errs() << space + "uint64_t max_RT = 0;\n";
        
        errs() << space + "for ( map<uint64_t, double>::reverse_iterator it = rt.rbegin(), eit = rt.rend(); it != eit; ++it) {\n";
        errs() << space + "    total_num_RT += it->second;\n";
        errs() << space + "    if (max_RT < it->first) {\n";
        errs() << space + "        max_RT = it->first;\n";
        errs() << space + "    }\n";
        errs() << space + "}\n";
        
        errs() << space + "double accumulate_num_RT = 0;\n";
        errs() << space + "for ( map<uint64_t, double>::reverse_iterator it = rt.rbegin(), eit = rt.rend(); it != eit; ++it) {\n";
        errs() << space + "    P[it->first] = accumulate_num_RT / total_num_RT;\n";
        errs() << space + "    accumulate_num_RT += it->second;\n";
        errs() << space + "}\n";
        
        errs() << space + "P[0] = 1;\n";
        
        errs() << space + "double sum_P = 0.0;\n";
        errs() << space + "uint64_t aet = 0, prev_aet = 0;\n";
        errs() << space + "while (sum_P < cache_size && aet <= max_RT) {\n";
        errs() << space + "    if (P.find(aet) != P.end()) {\n";
        errs() << space + "        sum_P += P[aet];\n";
        errs() << space + "        prev_aet = aet;\n";
        errs() << space + "    } else {\n";
        errs() << space + "        sum_P += P[prev_aet];\n";
        errs() << space + "    }\n";
        errs() << space + "    aet++;\n";
        errs() << space + "}\n";
        errs() << space + "return P[prev_aet];\n";
        errs() << "}\n";
        
        return;
    }

#if defined(REFERENCE_GROUP)
    void TuningCodeGen_ref::rtMergeGen() {
        errs() << "/* Merge the refRT to RT */\n";
        errs() << "void rtMerge() {\n";
        errs() << "    for(map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
        errs() << "        for (map<uint64_t, double>::iterator iit = it->second.begin(); iit != it->second.end(); ++iit) {\n"; 
        errs() << "            rtHistoCal(RT, iit->first, iit->second);\n";
        // errs() << "            subBlkRT(RT, iit->first, iit->second);\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        return;
    }
#endif

    void TuningCodeGen_ref::statDumpGen() {
        errs() << "/* Dump the reuse statistics */\n";
        errs() << "void statDump() {\n";
        errs() << "    cout << \"Total Neighboring Reuses: \" << total_neighbor / total_reuse << endl;\n";
        errs() << "    cout << \"Total Scaling Reuses: \" << total_scale / total_reuse << endl;\n";
        errs() << "    cout << \"Total Folding Src-Sink Reuses: \" << total_fold / total_reuse << endl;\n";
        errs() << "    cout << \"Total Inter Chunk Reuses: \" << total_interchunk / total_reuse << endl;\n";
        errs() << "    cout << \"Total Shared Reuses: \" << share_reuse / total_reuse << endl;\n";
        errs() << "    cout << \"Total Reuses: \" << total_reuse << endl;\n";
        errs() << "    return;\n";
        errs() << "}\n";
    }

    string TuningCodeGen_ref::getBound(Value *bound) {
        
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
    
    string TuningCodeGen_ref::getBound_Start(Value *bound) {
        
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

    string TuningCodeGen_ref::getLoopInc(Value *inc) {
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

    uint64_t TuningCodeGen_ref::getOuterMostLoopIterationSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node) {
        return ((uint64_t) stoi(getBound((*node->LIS->LB)[0].second)) - (uint64_t) stoi(getBound((*node->LIS->LB)[0].first))) / (uint64_t) stoi(getLoopInc((*node->LIS->INC)[0]));
    }

    uint64_t TuningCodeGen_ref::computeIterSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, bool inclusive) {
        uint64_t iterSpace = 0;
        if (LoopRefTree->next != NULL) {
            for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                if ((*it)->L == NULL) {
                    iterSpace += 1;
                } else {
                    iterSpace += computeIterSpace(*it, true);
                }
            } 
            if (inclusive) {
                iterSpace *= getOuterMostLoopIterationSpace(LoopRefTree);
            }
        }
        return iterSpace;
    }
    
    vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> TuningCodeGen_ref::findLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, string refName, int useID, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops) {
        
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
            for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopTmp = findLoops(*it, refName, useID, loops);
                if (loopTmp.size() != 0) {
                    loopRes = loopTmp;
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
    void TuningCodeGen_ref::searchReuseDifferentLoopsUpdateFuncGen() {
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
    
    void TuningCodeGen_ref::searchReuseDifferentLoopsCalFuncGen() {
        
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
    
    bool TuningCodeGen_ref::searchReuseDifferentLoopsInitGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, bool GenFlag,  string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> currentLoops, string space) {
        
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
            for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops_New = currentLoops;
                if (LoopRefTree->L != NULL) {
                    currentLoops_New.push_back(LoopRefTree);
                }
                if ((*it)->isThreadNode) { continue; }
                GenFlag = searchReuseDifferentLoopsInitGen(*it, GenFlag, refName, useID, loops, currentLoops_New, space);
            }
        }
        
        return GenFlag;
    }
    
    bool TuningCodeGen_ref::searchReuseDifferentLoopsBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag,  string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> currentLoops, string space) {
        
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
            for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
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
    void TuningCodeGen_ref::searchReuseSameLoopInitGen( string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, string space) {
        
        //errs() << "Search result reuse: (ref name " << refName << " ) ( ID " << useID << " ) ( numOfLoops " << loops.size() << " )\n";

        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator ait = loops.back()->next->begin(), eit = loops.back()->next->end(); ait != eit; ++ait) {
            if ((*ait)->isThreadNode) { continue; }
            if ((*ait)->AA != NULL) {
                if (arrayName[(*ait)->AA] == refName) {
                    errs() << space + "uint64_t prev_cnt_" + refName +  to_string(refNumber[(*ait)->AA]) + " = -1;\n";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        errs() << space + "uint64_t prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_Start_" + refName +  to_string(refNumber[(*ait)->AA]) + " = -1;\n";
                        errs() << space + "uint64_t prev_" + indvName[(*(*it)->LIS->IDV)[0]] + "_End_" + refName +  to_string(refNumber[(*ait)->AA]) + " = -1;\n";
                    }
                }
            }
        }
    
        return;
    }
    
    void TuningCodeGen_ref::searchReuseSameLoopBodyGen( string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops,  string space) {
        
        for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator ait = loops.back()->next->begin(), eit = loops.back()->next->end(); ait != eit; ++ait) {
            if ((*ait)->isThreadNode) { continue; }
            if ((*ait)->AA != NULL) {
                if (arrayName[(*ait)->AA] == refName) {
                    errs() << space + "if ( prev_cnt_" + refName +  to_string(refNumber[(*ait)->AA]) + " != -1) {\n";
                    errs() << space + "    if ( calAddr" + refName +  to_string(refNumber[(*ait)->AA]) + "( ";
                    string tmp = "";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            
                            tmp += indvName[(*(*it)->LIS->IDV)[i]] + "_Start";
                            tmp += " - prev_" + indvName[(*(*it)->LIS->IDV)[i]] + "_Start_" + refName +  to_string(refNumber[(*ait)->AA]);
                            tmp += " + prev_" + indvName[(*(*it)->LIS->IDV)[i]] + "_End_" + refName +  to_string(refNumber[(*ait)->AA]);
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
                    errs() << space + "        pred = " + "prev_cnt_" + refName +  to_string(refNumber[(*ait)->AA]) + ";\n";
#endif
                    
                    errs() << space + "        rtHistoCal(prev_cnt_" + refName +  to_string(refNumber[(*ait)->AA]);
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
    bool TuningCodeGen_ref::refRTSearchGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag, string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, string space) {
        
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
                            errs() << indvName[(*(*it)->LIS->IDV)[0]];
                            errs() << " == " + indvName[(*(*it)->LIS->IDV)[0]] + "_Start";
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
                    
                    errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] << "=" << getBound((*LoopRefTree->LIS->INC)[0]);
                    
                } else if ((*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGE || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGT || (*LoopRefTree->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGT) {
                    
                    errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] << "="<< getBound((*LoopRefTree->LIS->INC)[0]);
                    
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

                    string is_normal_ref = "true";
                    string is_in_same_loop = "true";
                    string outermost_src_indv = "(" + indvName[(*(*loops.begin())->LIS->IDV)[0]] + "_Start -" + getBound_Start((*(*loops.begin())->LIS->LB)[0].first) + ")";
                    string outermost_sink_indv = "(" + indvName[(*(*currentLoops.begin())->LIS->IDV)[0]] + " - " + getBound_Start((*(*currentLoops.begin())->LIS->LB)[0].first) + ")";
                    if (arrayTypeMap[refName + std::to_string(useID)] != REGULAR || arrayTypeMap[refName + std::to_string(refNumber[LoopRefTree->AA])] != REGULAR) {
                        is_normal_ref = "false";
                    } 
                    
                    bool is_src_loop_outermost = find(outloops.begin(), outloops.end(), (*loops.begin())) != outloops.end();
                    bool is_sink_loop_outermost = find(outloops.begin(), outloops.end(), (*currentLoops.begin())) != outloops.end();
                    errs() << space << "            /* is_src_loop_outermost: " << is_src_loop_outermost << " */\n";
                    errs() << space << "            /* is_sink_loop_outermost: " << is_sink_loop_outermost << " */\n";
                    /* source and sink are not in the same loop 
                     * 
                     */
                    if ((*currentLoops.begin()) != (*loops.begin())) {
                        is_in_same_loop = "false";
                    }
                    errs() << space << "            /* is_normal_ref: " << is_normal_ref << " */\n";
                    errs() << space << "            /* is_in_same_loop: " << is_in_same_loop << " */\n";
                    errs() << space << "            /* register the src/sink addr calculation function. Calling this function in other function do not have to know all index variables */\n";
                    errs() << space << "            function<uint64_t(uint64_t)> srcAddrCal = bind(calAddr" << refName <<  to_string(useID) << ", _1";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        if (find(outloops.begin(), outloops.end(), *it) != outloops.end()) { continue; }
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            errs() << ", ";
                            errs() << indvName[(*(*it)->LIS->IDV)[i]] << "_Start";
                        }
                    }
                    errs() << ");\n";
                    errs() << space << "            function<uint64_t(uint64_t)> sinkAddrCal = bind(calAddr" << refName << to_string(refNumber[LoopRefTree->AA]) << ", _1";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        if (find(outloops.begin(), outloops.end(), *it) != outloops.end()) { continue; }
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            errs() << ", ";
                            errs() << indvName[(*(*it)->LIS->IDV)[i]];
                        }
                    }
                    errs() << ");\n";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        errs() << space + "            /* ";
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            errs() << indvName[(*(*it)->LIS->IDV)[i]];
                            errs() << " ";
                        }
                        errs() << computeIterSpace(*it, false) << " */\n";
                    }
                    errs() << space + "            /* compute the number of accesses between source and sink chunk */\n";
                    /* first we compute middle chunks between source and sink chunk */
                    vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator head = find(outloops.begin(), outloops.end(), (*loops.begin()));
                    vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator end = find(outloops.begin(), outloops.end(), (*currentLoops.begin()));
                    uint64_t middle_access = 0;
                    errs() << space + "            uint64_t middle_accesses = 0;\n";
                    if (is_in_same_loop == "false") {
                        if (head != outloops.end() && end != outloops.end()) {
                            for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = head + 1; it < end; ++it) {
                                middle_access += computeIterSpace(*it, true);
                            }
                        }
                        errs() << space + "            if (isCompleteChunk(" << getOuterMostLoopIterationSpace((*loops.begin())) << ", thread_cnt)) {\n";
                        errs() << space << space << "            middle_accesses += ";
                        /* Here we compute the remaining chunks to be iterated after source chunk and chunks to be iterated before sink chunk. Then add the iterations computed in the first step */
                        errs() << "((";
                        if ((*currentLoops.begin()) != (*loops.begin())) {
                            errs() << "getChunkNum(" << getOuterMostLoopIterationSpace((*loops.begin())) << ", thread_cnt)";
                        }
                        errs() << " - getChunkID(" << outermost_src_indv << ", thread_cnt) - 1) * (CHUNK_SIZE * thread_cnt)) * " << to_string(outMostLoopPerIterationSpace[refName +  to_string(useID)]) << " + getChunkID(" << outermost_sink_indv << ", thread_cnt) * CHUNK_SIZE * thread_cnt * " << to_string(outMostLoopPerIterationSpace[refName +  to_string(refNumber[LoopRefTree->AA])]) << ";\n";
                        errs() << space << "            } else {\n";
                        errs() << space << space << "            middle_accesses += ";
                        /* Here we compute the remaining chunks to be iterated after source chunk and chunks to be iterated before sink chunk. Then add the iterations computed in the first step */
                        errs() << "((";
                        if ((*currentLoops.begin()) != (*loops.begin())) {
                            errs() << "getChunkNum(" << getOuterMostLoopIterationSpace((*loops.begin())) << ", thread_cnt)";
                        }
                        // -2 means remove the source and the last chunk
                        errs() << " - getChunkID(" << outermost_src_indv << ", thread_cnt) - 2) * CHUNK_SIZE * thread_cnt + (" << getOuterMostLoopIterationSpace((*loops.begin())) << " % (CHUNK_SIZE * thread_cnt))" << ") * " << to_string(outMostLoopPerIterationSpace[refName +  to_string(useID)]) << " + getChunkID(" << outermost_sink_indv << ", thread_cnt) * CHUNK_SIZE * thread_cnt * " << to_string(outMostLoopPerIterationSpace[refName +  to_string(refNumber[LoopRefTree->AA])]) << ";\n";
                        errs() << space << "            }\n";
                    }
                    errs() << space + "            middle_accesses += "<< middle_access << ";\n";
                    errs() << "#ifdef DEBUG\n";
                    errs() << space + "            cout << \" middle_access is \" << middle_accesses << endl;\n";
                    errs() << "#endif\n";
                    errs() << space + "            int reuse_type = -1;\n";
#ifdef SMOOTHING
                    errs() << space + "            uint64_t step = cnt;\n";
                    vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator src_it = loops.begin();
                    vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator sink_it = currentLoops.begin();
                    while(src_it != loops.end() && sink_it != currentLoops.end()) {
                        errs() << space + "            if ( " << indvName[(*(*src_it)->LIS->IDV)[0]] << "_Start != " << indvName[(*(*sink_it)->LIS->IDV)[0]] << ") {\n";
                        // errs() << space << "                step = " << computeIterSpace(*src_it, false) << ";\n";
                        if (find(outloops.begin(), outloops.end(), *sink_it) != outloops.end()) { 
                            errs() << space << "                step = " << computeIterSpace(*src_it, false) << ";\n";
                        } else {
                            errs() << space << "                step = abs(" << indvName[(*(*sink_it)->LIS->IDV)[0]] << " -" << indvName[(*(*src_it)->LIS->IDV)[0]] << "_Start" << ") *" << computeIterSpace(*src_it, false) << ";\n";
                        }
                        errs() << space + "            }\n";
                        src_it++;
                        sink_it++;
                    }
                    errs() << space + "            uint64_t parallel_rt = parallel_predict(" << outermost_src_indv << ", " << outermost_sink_indv << ", cnt, " << to_string(outMostLoopPerIterationSpace[refName +  to_string(useID)]) << ", " << to_string(outMostLoopPerIterationSpace[refName +  to_string(refNumber[LoopRefTree->AA])]) << ", middle_accesses, step, thread_cnt, " << is_normal_ref << ", " << is_in_same_loop << ", reuse_type, srcAddrCal, sinkAddrCal);\n";
#else
                    errs() << space + "            uint64_t parallel_rt = parallel_predict(" << outermost_src_indv << ", " << outermost_sink_indv << ", cnt, " << to_string(outMostLoopPerIterationSpace[refName +  to_string(useID)]) << ", " << to_string(outMostLoopPerIterationSpace[refName +  to_string(refNumber[LoopRefTree->AA])]) << ", middle_accesses, thread_cnt, " << is_normal_ref << ", " << is_in_same_loop << ", reuse_type, srcAddrCal, sinkAddrCal);\n";
#endif
#ifdef REFERENCE_GROUP
                    // errs() << space + "            refSubBlkRT(refRT, cnt, 1.0, \"" + refName + std::to_string(useID) + "\");\n";
                    errs() << space + "            refRTHistoCal(refRT, parallel_rt, 1.0, \"" + refName + std::to_string(useID) + "\"";
#else
                    // errs() << space + "            subBlkRT(RT, cnt, 1.0);\n";
                    errs() << space + "            if (parallel_rt == 0) { goto EndSample; }\n";
                    errs() << "#ifdef DEBUG\n";
                    errs() << space + "            rtHistoCal(RT, parallel_rt, 1.0);\n";
#endif
                    errs() << "#else\n";
                    errs() << space + "            rtHistoCal(RT, parallel_rt, 1.0);\n";
                    // errs() << space + "        subBlkRT(RT, parallel_rt, 1.0);\n";
                    errs() << "#endif\n";
                    errs() << "#ifdef DEBUG\n";
                    // errs() << space + "            if (parallel_rt != 1 && parallel_rt != 13) {\n";
                    errs() << space << "                cout << \"[" << refName << to_string(useID) << " --> " << refName << to_string(refNumber[LoopRefTree->AA]) << "]" << " [\" << parallel_rt << \"] (\" << ";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            errs() << indvName[(*(*it)->LIS->IDV)[i]] << "_Start";
                            if (it != loops.end()-1) {
                                errs() << "<< \", \" << ";
                            }
                        }
                    }
                    errs() << "<< \") --> (\" << ";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            errs() << indvName[(*(*it)->LIS->IDV)[i]];
                            if (it != currentLoops.end()-1) {
                                errs() << "<< \", \" << ";
                            }
                        }
                    }
                    errs() << "<< \") \" << endl;\n"; 
                    // errs() << space + "            }\n";
                    errs() << "#endif\n";
                    /* output format:
                    * src-ref, src-idx, sink-ref, sink-idx, ris, rip, type
                    */
                    errs() << space << "                // cout << \"" << refName << to_string(useID) << ";\" << ";
                    errs() << "\"(\" << ";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            errs() << indvName[(*(*it)->LIS->IDV)[i]] << "_Start";
                            if (it != loops.end()-1) {
                                errs() << "<< \", \" << ";
                            }
                        }
                    }
                    errs() << "<< \");\" << \"" << refName << to_string(refNumber[LoopRefTree->AA]) << ";\" << ";
                    errs() << "\"(\" << ";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            errs() << indvName[(*(*it)->LIS->IDV)[i]];
                            if (it != currentLoops.end()-1) {
                                errs() << "<< \", \" << ";
                            }
                        }
                    }
                    errs() << "<< \");\" << cnt << \";\" << parallel_rt << \";\" << reuse_type << endl;\n"; 
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
            for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops_New = currentLoops;
                if (LoopRefTree->L != NULL) {
                    currentLoops_New.push_back(LoopRefTree);
                }
                GenFlag = refRTSearchGen(*it, GenFlag, refName, useID, loops, currentLoops_New, space + "    ");
            }
            if (LoopRefTree->L != NULL && GenFlag == true) {
                errs() << space + "}\n";
                errs() << space + "}\n";
            }
        }
        
        return GenFlag;
    }


    void TuningCodeGen_ref::refRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, string refName, int useID) {
        
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

                errs() << space << "    if (" << indvName[(*(*lit)->LIS->IDV)[i]] << "_Start % " <<  getLoopInc((*(*lit)->LIS->INC)[i]) << " != 0) goto SAMPLE; \n";
                if (find(outloops.begin(), outloops.end(), *lit) != outloops.end()) {
                    errs() << space << "    if (" << indvName[(*(*lit)->LIS->IDV)[i]] << "_Start + thread_cnt * CHUNK_SIZE > " << getBound_Start((*(*lit)->LIS->LB)[i].second) << ") { goto SAMPLE; }\n";
                }
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
    
    void TuningCodeGen_ref::refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
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
            errs() << "void ref_" + arrayName[it->first] +  to_string(it->second) + "(int thread_cnt, map<uint64_t, double> & RT) {\n";
            refRTBodyGen(LoopRefTree, arrayName[it->first], it->second);
            errs() << "}\n";
        }
        
        return;
    }


    void TuningCodeGen_ref::tuningFuncGen() {
        string space = "    ";
        errs() << "void mergePrivateRT(map<int, map<uint64_t, double>> prt) {\n";
        errs() << space << "for(map<int, map<uint64_t, double>>::iterator it = prt.begin(); it != prt.end(); ++it) {\n";
        errs() << space << space << "map<uint64_t, double> tmp = it->second;\n";
        errs() << space << space << "for(map<uint64_t, double>::iterator iit = tmp.begin(); iit != tmp.end(); ++iit) {\n";
        errs() << space << space << space << "rtHistoCal(RT, iit->first, iit->second);\n";
        errs() << space << space << "}\n";
        errs() << space << "}\n";
        errs() << space << "return;\n";
        errs() << "}\n";
        errs() << "void generate_per_thread_reuse(int thread_cnt) {\n";
        errs() << space << "map<uint64_t, double> tmpRT;\n";
        errs() << space << "double accumulate_num_RT = 0.0;\n";
        errs() << space << "double expected_rt = 0.0;\n";
        for ( map<Instruction*, int>::iterator it = refNumber.begin(), eit = refNumber.end(); it != eit; ++it) {
#ifdef PARALLEL_CXX_THREAD
            /*
                errs() << space + " thread t_"+ (*it).first + "_" + to_string(i)+ "(";
                errs() << "ref_" + (*it).first +  to_string(i) + ");\n";
             */
            errs() << space + "thread t_"+ arrayName[it->first] + "_" +  to_string(it->second) + "(";
            errs() << "ref_" + arrayName[it->first] +  to_string(it->second) + ", thread_cnt, tmpRT);\n";
#else
            /*
                errs() << space + "ref_" + (*it).first +  to_string(i) + "();\n";
             */
            errs() << space << "ref_" + arrayName[it->first] +  to_string(it->second) + "(thread_cnt, tmpRT);\n";
#endif
        /*
            }
         */
        }
        errs() << space << "per_tcount_RT[thread_cnt] = tmpRT;\n";
        errs() << space << "/* iterate each tid and dump its LLC miss ratio (13.75MB) */\n";
        errs() << space << "per_tcount_MR[thread_cnt] = RTtoMR_AET_C(tmpRT, 2285280);\n";
        errs() << "}\n";
        errs() << "/* Compute expected reuse for current thread number */\n";
        errs() << "void compute_expected_reuse(int thread_cnt) {\n";
        errs() << space << "double accumulate_num_RT = 0.0;\n";
        errs() << space << "double expected_rt = 0.0;\n";
        errs() << space << "map<uint64_t, double> tmpRT = per_tcount_RT[thread_cnt];\n";
        errs() << space << "for ( map<uint64_t, double>::iterator it = tmpRT.begin(), eit = tmpRT.end(); it != eit; ++it) {\n";
        errs() << space << "    accumulate_num_RT += it->second;\n";
        errs() << space << "}\n";
        errs() << space << "for ( map<uint64_t, double>::iterator it = tmpRT.begin(), eit = tmpRT.end(); it != eit; ++it) {\n";
        errs() << space << "    expected_rt += (it->first * it->second / accumulate_num_RT);\n";
        errs() << space << "}\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << "cout << thread_cnt << \", \" << expected_rt << endl;\n";
        errs() << "#endif\n";
        errs() << space << "per_threadcnt_expected_rt[thread_cnt] = expected_rt;\n";
        errs() << "}\n";
    }
    
    void TuningCodeGen_ref::mainGen() {
    
        errs() << "int main() {\n";
        
        string space = "    ";
        /*
        for ( map<string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < (*it).second; i++) {
         */
#ifdef DumpRTMR
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "    PAPI_timer_init();\n";
        errs() << "    PAPI_timer_start();\n";
        errs() << "#endif\n";
#endif
        errs() << space << "int tlb = THREAD_LB;\n";
        errs() << space << "int tub = THREAD_UB;\n";
        errs() << space << "/* metadata to derive best thread */\n";
        errs() << space << "double min_expected_rt = 0.0;\n";
        errs() << space << "int min_tnum = tlb;\n";
        errs() << space << "vector<thread> thread_vec;\n";
        errs() << space << "for (int t = tlb; t <= tub; t++) {\n";
        errs() << space << space << "/* Currently we consider all even number threads only */\n";
        errs() << space << space << "if ( ceil(log2(t)) != floor(log2(t)) ) { continue; }\n";
        errs() << space << space << "generate_per_thread_reuse(t);\n";
        // errs() << space << space << "if (t == THREAD_NUM) {\n";
        // errs() << space << space << space << "mergePrivateRT(privateRT);\n";
        // errs() << space << space << space << "RTtoMR_AET();\n";
        // errs() << space << space << space << "dumpMR();\n";
        // errs() << space << space << "}\n";
        // errs() << space << space << "/* iterate each tid and dump its private L2 miss ratio (1MB) */\n";
        // errs() << space << space << "for (int id = 0; id < t; id++) {\n";
        // errs() << space << space << space << "cout << RTtoMR_AET_C(privateRT[id], 16384) << \" | \";\n";
        // errs() << space << space << "}\n";
        // errs() << space << space << "cout << endl;";
        errs() << space << space << "thread_vec.push_back(thread(compute_expected_reuse, t));\n";
        errs() << space << "}\n";
        errs() << space << "/* compute expected rt for each thread count */\n";
        errs() << space << "for (vector<thread>::iterator threadit = thread_vec.begin(); threadit != thread_vec.end(); ++threadit) {\n";
        errs() << space << space << "(*threadit).join();\n";
        errs() << space << "}\n";
        errs() << space << "for (map<int, double>::iterator mit = per_threadcnt_expected_rt.begin(); mit != per_threadcnt_expected_rt.end(); ++mit) {\n";
        // errs() << "#ifdef DEBUG\n";
        errs() << space << space << "cout << mit->first << \", \" << mit->second << endl;\n";
        // errs() << "#endif\n";
        errs() << space << space << "if (mit == per_threadcnt_expected_rt.begin()) {\n;";
        errs() << space << space << space << "min_expected_rt = mit->second;\n";
        errs() << space << space << space << "continue;\n";
        errs() << space << space << "}\n";
        errs() << space << space << "if (min_expected_rt > mit->second) {\n";
        errs() << space << space << space << "min_tnum = mit->first;\n";
        errs() << space << space << space << "min_expected_rt = mit->second;\n";
        errs() << space << space << "}\n";
        errs() << space << "}\n";
        // // errs() << space << "for (map<int, double>::iterator mit = per_tcount_MR.begin(); mit != per_tcount_MR.end(); ++mit) {\n";
        // // errs() << space << space << "cout << mit->second << \" \";\n";
        // // errs() << space << "}\n";
        /*
            }
         */
    
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
        errs() << "    cout << min_tnum << endl;\n";
#ifdef DumpRTMR
// #ifdef REFERENCE_GROUP
//         errs() << "    rtMerge();\n";
// #endif
        // errs() << "    RTtoMR_AET();\n";
        errs() << "#ifdef PAPI_TIMER\n";
        errs() << "    PAPI_timer_end();\n";
        errs() << "    PAPI_timer_print();\n";
        errs() << "#endif\n";
// #ifdef REFERENCE_GROUP
//         errs() << "    refRTDump();\n";
// #endif
//         errs() << "    cout << \"Samples: \" << \"" << l1_dcache_access << "\" << endl;\n"; 
//         errs() << "    rtDump();\n";
//         errs() << "    dumpMR();\n";
        // errs() << "    statDump();\n";
#else
        // errs() << "#ifdef PAPI_TIMER\n";
        // errs() << "// Get ending timepoint\n"; 
        // errs() << "    auto stop = high_resolution_clock::now(); \n";
  
        // errs() << "    // Get duration. Substart timepoints to\n";
        // errs() << "    // get durarion. To cast it to proper unit\n"; 
        // errs() << "    // use duration cast method\n";
        // errs() << "    auto duration = duration_cast<microseconds>(stop - start);\n ";
        // errs() << "    cout << \"Time taken by SPS:  \" << duration.count() << endl; \n";
        // errs() << "#endif\n";
#endif
        errs() << "    return 0;\n";
        errs() << "}\n";
        
        return;
    }

    void TuningCodeGen_ref::defineRefs(Instruction* arrayInstr) {
        errs() << "/* Array " << arrayName[arrayInstr] << "\t";
        for(vector<string>::iterator vit = arrayAccessVariable[arrayInstr].begin(); vit != arrayAccessVariable[arrayInstr].end(); ++vit) {
            errs() << *vit << " ";
        }
        errs() << "*/ \n";
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* target = NULL;
        // find the outermost loop that this array enclosed in 
        for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = outloops.begin(); it != outloops.end(); ++it) {
            if (find(refPerOutMostLoop[*it].begin(), refPerOutMostLoop[*it].end(), arrayInstr) != refPerOutMostLoop[*it].end()) {
                target = *it;
                break;
            }
        }
        // this array expression does not in an loop
        if (target == NULL) {
            arrayTypeMap[arrayName[arrayInstr] + to_string(refNumber[arrayInstr])] = OUTERMOST_INVARIANT;
            return;
        }
        // for each induction variable in an array expression
        for(vector<string>::iterator vit = arrayAccessVariable[arrayInstr].begin(); vit != arrayAccessVariable[arrayInstr].end(); ++vit) {
            errs() << "/* " << indvName[(*(target)->LIS->IDV)[0]] << " */\n"; 
            // contains outermost loop induction variable
            if (indvName[(*(target)->LIS->IDV)[0]] == *vit ) {
                // contains only 1 induction 
                if (arrayAccessVariable[arrayInstr].size() == 1) {
                    arrayTypeMap[arrayName[arrayInstr] + to_string(refNumber[arrayInstr])] = OUTERMOST_ONLY;
                    return;
                } else if (vit != arrayAccessVariable[arrayInstr].begin()) {
                    /* in this case
                     * array index expression contains the outermost loop
                     * induction variable but it does not in the same access
                     * level as the outermost loop does.
                     * 
                     * i.e. i is the outermost loop and j is its inner loop
                     * A[j][i] meet this condition
                     */
                    errs() << "/* " << *vit << ", " << *(arrayAccessVariable[arrayInstr].begin()) << " */\n";
                    arrayTypeMap[arrayName[arrayInstr] + to_string(refNumber[arrayInstr])] = REVERSE;
                    return;
                }
                arrayTypeMap[arrayName[arrayInstr] + to_string(refNumber[arrayInstr])] = REGULAR;
                return;
            }
        }
        arrayTypeMap[arrayName[arrayInstr] + to_string(refNumber[arrayInstr])] = OUTERMOST_INVARIANT;
        return;
    }
    
    void TuningCodeGen_ref::parallelModelCodeGen() {
        string space = "    ";
        errs() << "bool isCompleteChunk(uint64_t is, int thread_cnt) {\n";
        errs() << space << "return ((is % (CHUNK_SIZE * thread_cnt)) == 0);\n";
        errs() << "}\n";

        errs() << "int getChunkNum(uint64_t is, int thread_cnt) {\n";
        errs() << space << "if (is % (CHUNK_SIZE * thread_cnt) != 0) {\n";
        errs() << space << space << "return is / (CHUNK_SIZE * thread_cnt) + 1;\n";
        errs() << space << "}\n";
        errs() << space << "return is / (CHUNK_SIZE * thread_cnt);\n";
        errs() << "}\n";

        errs() << "int getChunkID(uint64_t i, int thread_cnt) {\n";
        errs() << space << "return floor(i / (CHUNK_SIZE * thread_cnt));\n";
        errs() << "}\n";

        errs() << "int getThreadID(uint64_t i, int thread_cnt) {\n";
        errs() << space << "return i / CHUNK_SIZE - floor(i / (CHUNK_SIZE * thread_cnt))*thread_cnt ;\n";
        errs() << "}\n";

        errs() << "int getThreadLocalPos(uint64_t i) {\n";
        errs() << space << "return i \% CHUNK_SIZE;\n";
        errs() << "}\n";
        errs() << "int search_src_candidate_neighbor(uint64_t i, int thread_cnt, function<uint64_t(uint64_t)> calAddr) {\n";
        errs() << space << "int c_start = i \% CHUNK_SIZE + i / (CHUNK_SIZE * thread_cnt) * thread_cnt * CHUNK_SIZE;\n";
        errs() << space << "int c_end = c_start + (thread_cnt - 1) * CHUNK_SIZE;\n";
        errs() << space << "for (int c = i + CHUNK_SIZE; c <= c_end; c=c+CHUNK_SIZE) {\n";
        errs() << space << space << "if (calAddr(i) == calAddr(c)) { return getThreadID(c, thread_cnt); }\n";
        errs() << space << "}\n";
        errs() << space << "return -1;\n";
        errs() << "}\n";
        errs() << "int search_sink_candidate_neighbor(uint64_t i, int thread_cnt, function<uint64_t(uint64_t)> calAddr) {\n";
        errs() << space << "int c_start = i \% CHUNK_SIZE + i / (CHUNK_SIZE * thread_cnt) * thread_cnt * CHUNK_SIZE;\n";
        errs() << space << "for (int c = c_start; c <= i; c=c+CHUNK_SIZE) {\n";
        errs() << space << space << "if (calAddr(i) == calAddr(c)) { return getThreadID(c, thread_cnt); }\n";
        errs() << space << "}\n";
        errs() << space << "return -1;\n";
        errs() << "}\n";
#ifdef SMOOTHING
        errs() << "void parallel_smoothing(uint64_t ris, uint64_t rip,  uint64_t l, int thread_cnt, ReuseType type) {\n";
        errs() << space << "if (type == ReuseType::INTER_CHUNK) {\n";
        errs() << space << space << "total_interchunk += 1.0;\n";
        errs() << space << "}\n";
        errs() << space << "else if (type == ReuseType::NEIGHBOR) {\n";
        errs() << space << space << "vector<double> constant = { 0.51, 0.38, 0.10, 0.01 };\n";
        // errs() << space << space << "vector<double> constant = { 0.25, 0.25, 0.25, 0.25 };\n"; 
        // errs() << space << space << "if (ris < 64) { constant = {0.25, 0.25, 0.25, 0.25}; }\n";
        errs() << space << space << "if (ris == 1) {\n";
        errs() << space << space << space << "rtHistoCal(RT, ris, 1.0);\n";
        errs() << space << space << space << "return;\n";
        errs() << space << space << "}\n";
        // errs() << space << "if (ri >= 64) {\n";
        errs() << space << space << "/*\n";
        errs() << space << space << "for (int t = 0; t < thread_cnt; t++) {\n";
        // errs() << space << space << space << "uint64_t lb = rip + (t-1) * l >= 1 ? rip + (t-1) * l : 1;\n";
        // errs() << space << space << space << "for (uint64_t i = lb; i <= rip + (t * l); i++) {\n";
        errs() << space << space << space << "uint64_t lb = ris + (t-1) * l >= 1 ? ris + (t-1) * l : 1;\n";

        errs() << space << space << space << "for (uint64_t i = lb; i <= ris + (t * l); i++) {\n";
        errs() << space << space << space << space << "rtHistoCal(RT, i, constant[t] / l);\n";
        errs() << space << space << space << "}\n";
        errs() << space << space << "}\n";
        errs() << space << space << "*/\n";
        errs() << space << space << "for (int i = 0; i < constant.size(); i++) {\n";
        errs() << space << space << space << "rtHistoCal(RT, ris * pow(2, i-1), constant[i]);\n";
        errs() << space << space << "}\n";
        errs() << space << "}\n";
        errs() << space << "else if (type == ReuseType::SCALE) {\n";
        errs() << space << space << "vector<double> constant = { 0.13, 0.18, 0.44, 0.25 };\n";
        // errs() << space << space << "if (ris < 64) { constant = {0.25, 0.25, 0.25, 0.25}; }\n";
        errs() << space << space << "for (int i = 0; i < constant.size(); i++) {\n";
        errs() << space << space << space << "rtHistoCal(RT, ris * pow(2, i), constant[i]);\n";
        errs() << space << space << "}\n";
        errs() << space << space << "/*\n";
        errs() << space << space << "for (int t = 0; t < thread_cnt; t++) {\n";
        errs() << space << space << space << "for (uint64_t i = ris + (pow(2, t) * l); i < ris + (pow(2, (t+1)) * l); i++) {\n";
        errs() << space << space << space << space << "rtHistoCal(RT, i, constant[t] / l * (pow(2, t+1) - pow(2, t)));\n";
        errs() << space << space << space << "}\n";
        errs() << space << space << "}\n";
        errs() << space << space << "*/\n";
        errs() << space << "}\n";
        errs() << space << "else if (type == ReuseType::FOLD_SRC_SINK) {\n";
        errs() << space << space << "vector<double> constant = { 0.25, 0.25, 0.25, 0.25 };\n";
        // errs() << space << "if (ri >= 64) {\n";
        errs() << space << space << "for (int t = 0; t < thread_cnt; t++) {\n";
        errs() << space << space << space << "uint64_t lb = ris + (t * l) >= 1 ? ris + (t+1) * l : 1;\n";
        errs() << space << space << space << "for (uint64_t i = lb; i < (ris + t * l); i++) {\n";
        errs() << space << space << space << space << "rtHistoCal(RT, i, constant[t] / l);\n";
        errs() << space << space << space << "}\n";
        errs() << space << space << "}\n";
        errs() << space << "}\n";
        // errs() << space << "}\n";
        errs() << "}\n";
#endif
#ifdef SMOOTHING
        errs() << "uint64_t parallel_predict(uint64_t i_src, uint64_t i_sink, uint64_t rt, uint64_t lsrc, uint64_t lsink, uint64_t middle_accesses, uint64_t step, int thread_cnt, bool is_normal_ref, bool is_in_same_loop, function<uint64_t(uint64_t)> srcAddrCal, function<uint64_t(uint64_t)> sinkAddrCal) {\n";
#else
        errs() << "uint64_t parallel_predict(uint64_t i_src, uint64_t i_sink, uint64_t rt, uint64_t lsrc, uint64_t lsink, uint64_t middle_accesses, int thread_cnt, bool is_normal_ref, bool is_in_same_loop, int & type, function<uint64_t(uint64_t)> srcAddrCal, function<uint64_t(uint64_t)> sinkAddrCal) {\n";
#endif
        errs() << space << "total_reuse += 1.0;\n";
        errs() << space << "uint64_t parallel_rt = rt;\n";
        errs() << space << "int tsrc = getThreadID(i_src, thread_cnt);\n";
        errs() << space << "int tsink = getThreadID(i_sink, thread_cnt);\n";
        errs() << space << "int dT = tsink - tsrc;\n";
        errs() << space << "int sink_neighbor_delta = 0;\n";
        errs() << space << "if (!is_in_same_loop || getChunkID(i_src, thread_cnt) != getChunkID(i_sink, thread_cnt)) {\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << "cout << \"Inter Chunk Reuse\" << endl;\n";
        errs() << space << space << ";\n";
        errs() << space << space << "cout << \"rt \" << rt << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << "type = 0; // code for inter chunk reuse\n";
#ifdef SMOOTHING
        errs() << space << space << "parallel_smoothing(rt, rt * thread_cnt - CHUNK_SIZE * thread_cnt * (lsrc*(thread_cnt - tsrc) + lsink * tsink) + CHUNK_SIZE * thread_cnt * lsrc - (thread_cnt - 1) * middle_accesses + dT, step, ReuseType::FOLD_SRC_SINK);\n";
        errs() << space << space << space << "return 0;\n";
#else
        errs() << space << space << "parallel_rt = rt * thread_cnt - CHUNK_SIZE * thread_cnt * (lsrc*(thread_cnt - tsrc) + lsink * tsink) + CHUNK_SIZE * thread_cnt * lsrc - (thread_cnt - 1) * middle_accesses + dT;\n";
#endif
        errs() << space << space << "if (dT != 0) {\n";
        errs() << space << space << space << "share_reuse += 1.0;\n";
        errs() << space << space << "} else {\n";
        errs() << space << space << space << "uint64_t private_rt = (CHUNK_SIZE - getThreadLocalPos(i_src)) * lsrc + (getThreadLocalPos(i_sink) + 1) * lsink + (middle_accesses / thread_cnt );\n";
        errs() << space << space << space << "privateRTHistoCal(tsrc, privateRT, private_rt, 1.0);\n";
        errs() << space << space << space << "total_private_stat += 1.0;\n";
        errs() << space << space << "}\n";
        errs() << space << space << "return parallel_rt;\n";
        errs() << space << "} else if (!is_normal_ref) {\n";
        errs() << space << space << "/* intra chunk reuse */\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << "cout << \"Neighboring Effect\" << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << "if (tsrc == tsink) {\n";
        errs() << space << space << space << "privateRTHistoCal(tsrc, privateRT, rt, 1.0);\n";
        errs() << space << space << space << "total_private_stat += 1.0;\n";
        errs() << space << space << "}\n";
        errs() << space << space << "int tsrc_neighbor = search_src_candidate_neighbor(i_src, thread_cnt, srcAddrCal);\n";
        errs() << space << space << "int tsink_neighbor = search_sink_candidate_neighbor(i_sink, thread_cnt, sinkAddrCal);\n";
        errs() << space << space << "if (tsrc_neighbor >= 0) {\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Find sink in src neighbor at \" << tsrc_neighbor << endl;\n";
        errs() << space << space << space << "cout << \"Neighbor Effect: \" << tsrc_neighbor - tsrc << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "total_neighbor += 1.0;\n";
        errs() << space << space << space << "type = 1; // code for src neighboring effect\n";
        
#ifdef SMOOTHING
        errs() << space << space << space << "parallel_smoothing(rt, tsrc_neighbor - tsrc, step, ReuseType::NEIGHBOR);\n";
        errs() << space << space << space << "return 0;\n";
#else
        errs() << space << space << space << "return tsrc_neighbor - tsrc;\n";
#endif
        errs() << space << space << "} else if (tsink_neighbor >= 0 && tsink_neighbor != tsink) {\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Find sink in sink neighbor at \" << tsink_neighbor << endl;\n";
        errs() << space << space << space << "cout << \"Neighbor Effect: \" << rt * thread_cnt + tsink_neighbor - tsink << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "if (getChunkID(i_src, thread_cnt) == getChunkID(i_sink, thread_cnt)) {\n"; 
        
#ifdef SMOOTHING
        errs() << space << space << space << space << "parallel_smoothing(rt, rt * thread_cnt + tsink_neighbor - tsink, step, ReuseType::NEIGHBOR);\n";
        errs() << space << space << space << space << "return 0;\n";
#else
        errs() << space << space << space << space << "type = 2; // code for sink neighboring effect\n";
        errs() << space << space << space << space << "sink_neighbor_delta = tsink_neighbor - tsink;\n";
        // errs() << space << space << space << space << "return rt * thread_cnt + tsink_neighbor - tsink;\n";
#endif
        errs() << space << space << space << space << "total_neighbor += 1.0;\n";
        errs() << space << space << space << space << "share_reuse += 1.0;\n";
        errs() << space << space << space << "}\n";
        errs() << space << space << "}\n";
        errs() << space << "}\n";
        errs() << space << "if (getChunkID(i_src, thread_cnt) == getChunkID(i_sink, thread_cnt)) {\n"; 
        errs() << space << space << "/* same thread -- scaling effect */\n";
        errs() << space << space <<"if (dT == 0) {\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Scaling Effect\" << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "total_scale += 1.0;\n";
#ifdef SMOOTHING
        errs() << space << space << space << "parallel_smoothing(rt, rt * thread_cnt, rt, ReuseType::SCALE);\n";
        errs() << space << space << space << "return 0;\n";
#else
        errs() << space << space << space << "parallel_rt = rt * thread_cnt;\n";
        errs() << space << space << space << "type = 3; // code for scaling effect\n";
#endif
        errs() << space << space << space << "privateRTHistoCal(tsrc, privateRT, rt, 1.0);\n";
        errs() << space << space << space << "total_private_stat += 1.0;\n";
        errs() << space << space << "} else if (getThreadLocalPos(i_src) <= getThreadLocalPos(i_sink)) { // src-sink order\n";
        errs() << space << space << space << "if ((rt * thread_cnt - CHUNK_SIZE * lsrc * thread_cnt * dT + dT) < 0) { printf(\"NORMAL ORDER NEGATIVE PRI\\n\"); }\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Src-Sink Order Folding Effect\" << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "type = 4; // code for src-sink order folding effect\n";
        errs() << space << space << space << "total_fold += 1.0;\n";
#ifdef SMOOTHING
        errs() << space << space << space << "parallel_smoothing(rt, rt * thread_cnt - CHUNK_SIZE * lsrc * thread_cnt * dT + abs(dT), step, ReuseType::FOLD_SRC_SINK);\n";
        errs() << space << space << space << "return 0;\n";
#else
        errs() << space << space << space << "parallel_rt = rt * thread_cnt - CHUNK_SIZE * lsrc * thread_cnt * dT + abs(dT);\n";
#endif
        errs() << space << space << "} else { // sink-src order\n";
        errs() << space << space << space << "if ((rt * thread_cnt - CHUNK_SIZE * lsrc * thread_cnt * dT + dT) < 0) { printf(\"REVERSE ORDER NEGATIVE PRI\\n\"); }\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Sink-Src Order Folding Effect\" << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "// parallel_rt = CHUNK_SIZE * lsrc * thread_cnt * dT - (rt * thread_cnt) - abs(dT);\n";
        errs() << space << space << space << "type = 5; // code for sink-src order folding effect\n";
        errs() << space << space << space << "return 0;\n";
        errs() << space << space << "}\n";
        errs() << space << "}\n";
        errs() << space << "return parallel_rt + sink_neighbor_delta;\n";
        errs() << "}\n";
    }


    void TuningCodeGen_ref::SmoothingGen() {
        errs() << "/* Smoothing the per-reference RTHisto Uniformly. Equally splite the RT to a range of bins */\n";
        errs() << "void uniform_smoothing(map<uint64_t, double> &rth, int thread_cnt) {\n";
        errs() << "    map<uint64_t, double> tmp;\n";
        errs() << "    for(map<uint64_t, double>::iterator it = rth.begin(); it != rth.end(); ++it) {\n";
        errs() << "        uint64_t mu = it->first;\n";
        errs() << "        /* Do uniform distribution for all ri, distribute from ri to ri * thread_cnt  */\n";
        errs() << "        if (it->second <= 0.0) { continue; }\n";
        errs() << "        uint64_t start_b = mu;\n";
        errs() << "        uint64_t end_b = mu * thread_cnt;\n";
        errs() << "        double split_val = it->second / (end_b - start_b);\n";
        errs() << "        for (int b = start_b; b <= end_b; b++) {\n";
        errs() << "            rtHistoCal(tmp, b, split_val);\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    rth.clear();\n";
        errs() << "    for(map<uint64_t, double>::iterator it = tmp.begin(); it != tmp.end(); ++it) {\n";
        // errs() << "        subBlkRT(rth, it->first, it->second);\n";
        errs() << "        rtHistoCal(rth, it->first, it->second);\n";
        // errs() << "        RT[it->first] = it->second;\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        return;
    }
    
    bool TuningCodeGen_ref::runOnFunction(Function &F) {
    
        errs() << " // Start to generating Static Sampling Code (reference based)\n";
        
        /* reading info from previous passes */
        arrayName = getAnalysis<idxAnalysis::IndexAnalysis>().arrayName;
        arrayExpression = getAnalysis<idxAnalysis::IndexAnalysis>().arrayExpression;
        arrayAccessVariable= getAnalysis<idxAnalysis::IndexAnalysis>().arrayAccessVariable;
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree;
        sampleNum = getAnalysis<sampleNumAnalysis::SampleNumberAnalysis>().sampleNum;
        l1_dcache_access = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().total_cache_access;
#ifdef PARALLEL
        map<Instruction*, uint64_t> tempL = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().outMostLoopPerIterationSpace;;
        outloops = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().outMostLoops;
        refPerOutMostLoop = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().refPerOutMostLoop;
#endif

        /* init */
        initArrayName();
        numberRefToSameArray(LoopRefTree);
        numberLoops(LoopRefTree);
        initIndvName(LoopRefTree);
#ifdef PARALLEL
        for(map<Instruction*, uint64_t>::iterator it = tempL.begin(); it != tempL.end(); ++it) {
            errs() << "/* " << arrayName[(it->first)] << to_string(refNumber[(it->first)]) << "\t" << it->second << " */\n";
            outMostLoopPerIterationSpace[arrayName[(it->first)] + to_string(refNumber[(it->first)])] = it->second;
        } 
#endif
        /* generate headers */
        headerGen();

        /* generate rtHistoCal function */
        rtHistoGen();

        /* generate subBlkRT function */
        subBlkRTGen();

        /* generate parallel predict model */
        parallelModelCodeGen();

#if defined(REFERENCE_GROUP)
        rtMergeGen();
#endif
        
        /* generate rtToMR function */
        rtToMRGen();
        
        /* generate rtDump function */
        // rtDumpGen();
        
        /* generate mrDump function */
        mrDumpGen();

        /* generate statDump function */
        statDumpGen();

#if defined(SMOOTHING)
        SmoothingGen();
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

        tuningFuncGen();
        
        /* generate main function */
        mainGen();

        return false;
    }
    
    void TuningCodeGen_ref::getAnalysisUsage(AnalysisUsage &AU) const {
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
