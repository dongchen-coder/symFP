#include "riCodeGen_ref.hpp"

namespace riCodeGen_ref {
    char AllLevelRICodeGen_ref::ID = 0;
    static RegisterPass<AllLevelRICodeGen_ref> X("riCodeGen_ref", "parallel random interleaving static sampling code generating pass (both reference and iteration)", false, false);

    AllLevelRICodeGen_ref::AllLevelRICodeGen_ref() : FunctionPass(ID) {}

    void AllLevelRICodeGen_ref::numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
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
    
    void AllLevelRICodeGen_ref::numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree->L != NULL) {
            if (loopNumber.find(LoopRefTree->L) == loopNumber.end()) {
                loopNumber[LoopRefTree->L] = loopNumber.size();
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                numberLoops(*it);
            }
        }
        return;
    }
    
    void AllLevelRICodeGen_ref::initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree == NULL) {
            return;
        }
        
        if (LoopRefTree->L != NULL) {
            for (std::vector<Value *>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvName[*it] = (*it)->getName();
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
    
    void AllLevelRICodeGen_ref::initArrayName() {
        
        for (std::map<Instruction*, std::string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            it->second.replace(std::find(it->second.begin(), it->second.end(), '.'), std::find(it->second.begin(), it->second.end(), '.') +1, 1, '_');
        }
        
        return;
    }
    
    void AllLevelRICodeGen_ref::addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<string> indvs) {
        
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
    
    void AllLevelRICodeGen_ref::addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        std::vector<string> indvs;
        addrCalFuncGen(LoopRefTree, indvs);
        
        return;
    }
    
    void AllLevelRICodeGen_ref::headerGen() {
        
        errs() << "#include <map>\n";
        errs() << "#include <set>\n";
        errs() << "#include <vector>\n";
        errs() << "#include <tuple>\n";
        errs() << "#include <algorithm>\n";
        errs() << "#include <cstdlib>\n";
        errs() << "#include <iostream>\n";
        errs() << "#include <cmath>\n";
#ifdef PARALLEL_CXX_THREAD
        errs() << "#include <thread>\n";
        errs() << "#include <mutex>\n";
#endif
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
        //errs() << "std::map<uint64_t, tuple<uint64_t, int>> LAT;\n";
        vector<string> visit;
        for (std::map<Instruction*, int>::iterator it = refNumber.begin(), eit = refNumber.end(); it != eit; ++it) {
            if (find(visit.begin(), visit.end(), arrayName[it->first]) != visit.end()) { continue; }
            visit.push_back(arrayName[it->first]);
            errs() << "std::map<uint64_t, tuple<uint64_t, int>> LAT_" + arrayName[it->first] + ";\n";
        }
        
        errs() << "std::map<uint64_t, double> RT;\n";
        errs() << "std::map<uint64_t, double> MR;\n";
        return;
    }

#ifdef DumpRTMR
    void AllLevelRICodeGen_ref::rtHistoGen() {
        errs() << "void rtHistoCal( map<uint64_t, double> &rth, int rt, int val ) {\n";
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
        errs() << "    if (rth.find(rt) == rth.end()) { \n";
        errs() << "        rth[rt] = 1;\n";
        errs() << "    } else {\n";
        errs() << "        rth[rt] += 1;\n";
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
        errs() << "void refRTHistoCal( int rt, int val, string ref ) {\n";
        errs() << "    if ( val <= 0) {\n;";
        errs() << "        return;\n";
        errs() << "    }\n";
        errs() << "    if (refRT.find(ref) == refRT.end()) { \n";
        errs() << "        refRT[ref][rt] = val;\n";
        errs() << "    }\n";
        errs() << "    else {\n";
        errs() << "        if (refRT[ref].find(rt) == refRT[ref].end()) { \n";
        errs() << "            refRT[ref][rt] = val;\n";
        errs() << "        } else {\n";
        errs() << "            refRT[ref][rt] += val;\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
#endif
    }
#endif

    /* Generate the function to calculate the bins */
    /* Generate the function to calculate the bins */
    void AllLevelRICodeGen_ref::subBlkRTGen() {
        string space = "    ";
        errs() << "void subBlkRT(map<uint64_t, double> &rth, int rt) {\n";
#ifdef PARALLEL_CXX_THREAD
        errs() << space + "unique_lock< mutex> lck (mtx, defer_lock);\n";
        errs() << space + "lck.lock();\n";
#endif
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
        errs() << space + "            rtHistoCal(rth, b - diff, 1);\n";
        errs() << space + "            break;\n";
        errs() << space + "        }\n";
        errs() << space + "    }\n";
        errs() << space + "}\n";
        errs() << space + "else {\n";
        errs() << space + "    rtHistoCal(rth, pow(2, msb-1), 1);\n";
        errs() << space + "}\n";
#ifdef PARALLEL_CXX_THREAD
        errs() << space + "lck.unlock();\n";
#endif
        errs() << space + "return;\n";
        errs() << "}\n";
#ifdef REFERENCE_GROUP  
        errs() << "\n";
        errs() << "void refSubBlkRT(int rt, string ref) {\n";
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
        errs() << space + "            refRTHistoCal(b - diff, 1, ref);\n";
        errs() << space + "            break;\n";
        errs() << space + "        }\n";
        errs() << space + "    }\n";
        errs() << space + "}\n";
        errs() << space + "else {\n";
        errs() << space + "    refRTHistoCal(pow(2, msb-1), 1, ref);\n";
        errs() << space + "}\n";
        errs() << space + "return;\n";
        errs() << "}\n";
#endif
        return;
    }
    
    
#ifdef DumpRTMR
    void AllLevelRICodeGen_ref::rtDumpGen() {
        
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
    
    void AllLevelRICodeGen_ref::rtToMRGen() {
        
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
    
    void AllLevelRICodeGen_ref::mrDumpGen() {
        
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
#endif

#if defined(PARALLEL) && defined(REFERENCE_GROUP)
    void AllLevelRICodeGen_ref::rtMergeGen() {
        errs() << "/* Merge the refRT to RT */\n";
        errs() << "void rtMerge() {\n";
        errs() << "    for(map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
        errs() << "        for (map<uint64_t, double>::iterator iit = it->second.begin(); iit != it->second.end(); ++iit) {\n"; 
        errs() << "            rtHistoCal(iit->first, iit->second);\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        return;
    }

    void AllLevelRICodeGen_ref::GaussianDistrGen() {
        errs() << "/* Distribute the per-reference RTHisto based on Gasussian */\n";
        errs() << "void gaussian_distr() {\n";
        errs() << "    for(map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
        errs() << "        map<uint64_t, double> tmp = it->second;\n";
        errs() << "        for (map<uint64_t, double>::iterator iit = tmp.begin(); iit != tmp.end(); ++iit) {\n"; 
        errs() << "            uint64_t mu = iit->first;\n";
        errs() << "            uint64_t Y = iit->second;\n";
        errs() << "            double sigma = 3 * (double)mu / 7;\n";
        errs() << "            uint64_t start_b = log2(mu) <= log2(THREAD_NUM) ? 0 : (log2(mu) - log2(THREAD_NUM));\n";
        errs() << "            uint64_t end_b = log2(mu) + log2(THREAD_NUM);\n";
        errs() << "            /* Clear the original refRT[ref]. Make sure that this will be done only once */\n";
        errs() << "            if (iit == tmp.begin()) {\n";
        errs() << "                refRT[it->first].clear();\n";
        errs() << "            }\n";
        errs() << "            for(int b = start_b; b <= end_b; b++) {\n";
        errs() << "                if (b >= 0) {\n";
        errs() << "                    double val = 0.0;\n";
        errs() << "                    val = (1 /  (sqrt(2 * M_PI) * sigma))* exp( -1 * pow((pow(2.0,(double)b) - mu), 2.0) / (2 * pow(sigma, 2.0)));\n";
        errs() << "                    cout << pow(2.0,(double)b) << \", \" << (val * Y) << endl;\n";
        errs() << "                    refRTHistoCal(pow(2.0,(double)b), (val * Y), it->first);\n";
        errs() << "                }\n";
        errs() << "            }\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        return;
        
    }

    void AllLevelRICodeGen_ref::UniformDistrGen() {
        errs() << "/* Distribute the per-reference RTHisto Uniformly. Equally splite the RT to a range of bins */\n";
        errs() << "void uniform_distr() {\n";
        errs() << "    for(map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
        errs() << "        map<uint64_t, double> tmp = it->second;\n";
        errs() << "        for (map<uint64_t, double>::iterator iit = tmp.begin(); iit != tmp.end(); ++iit) {\n";
        errs() << "            uint64_t mu = iit->first;\n";
        errs() << "            uint64_t start_b = log2(mu) <= log2(THREAD_NUM) ? 0 : (log2(mu) - log2(THREAD_NUM));\n";
        errs() << "            uint64_t end_b = (log2(mu) + log2(THREAD_NUM));\n";
        errs() << "            uint64_t split_val = iit->second / (end_b - start_b + 1);\n";
        errs() << "            /* Clear the original refRT[ref]. Make sure that this will be done only once */\n";
        errs() << "            if (iit == tmp.begin()) {\n";
        errs() << "                refRT[it->first].clear();\n";
        errs() << "            }\n";
        errs() << "            for(int b = start_b; b <= end_b; b++) {\n";
        errs() << "                if (b >= 0) {\n";
        errs() << "                    refRTHistoCal(pow(2.0, b), split_val, it->first);\n";
        errs() << "                }\n";
        errs() << "            }\n";
        errs() << "        }\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        return;
    }
#endif

    string AllLevelRICodeGen_ref::getBound(Value *bound) {
        
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
    
    string AllLevelRICodeGen_ref::getBound_Start(Value *bound) {
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
    
    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> AllLevelRICodeGen_ref::findLoops(
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, 
        //string refName,
        //int useID,
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops,
        bool enableOPT
    ) {
        
        /* If this is a loop, add it to the list */
        if (LoopRefTree->L != NULL) {
            loops.push_back(LoopRefTree);
        }
        /* If this is a reference and it's the reference we'd like to sample, return the list */
        if (LoopRefTree->AA != NULL) {
            /*
             if (arrayName[LoopRefTree->AA] == refName) {
                if (!enableOPT && refNumber[LoopRefTree->AA] == useID) {
                    return loops;
                } else if (enableOPT && refNumber[LoopRefTree->AA] >= useID) {
                    return loops;
                }
            } else {
                return std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>();
            }*/
            return loops;
        }
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopRes;
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                //std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopTmp = findLoops(*it, refName, useID, loops, enableOPT);
                std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopTmp = findLoops(*it, loops, enableOPT);
                if (loopTmp.size() != 0) {
                    loopRes = loopTmp;
                }
            }
        }
        return loopRes;
    }

    /* Generate the loop incremental func */
    void AllLevelRICodeGen_ref::LoopIterIncGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node,
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops,
            string space,
            bool allLevelInc
        ) {
        // getBound_Start((*(*lit)->LIS->LB)[0].first)
        errs() << space + "/* Iteration incrementation " << currentLoops.size() << " */\n";
        errs() << space + "progress[t_select][" + to_string(currentLoops.size()) + "] = progress[t_select][" + to_string(currentLoops.size()) + "] + 1; \n";
        if (allLevelInc) {
            int i = currentLoops.size()-1;
            for(vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::reverse_iterator rit = currentLoops.rbegin(); rit != currentLoops.rend(); ++rit) {
                if (i == currentLoops.size()-1) {
                    /* increment the index in the middle level */
                    errs() << space + "progress[t_select][" + to_string(i) + "] = progress[t_select][" + to_string(i) + "] + (progress[t_select][" + to_string(i+1) + "] / " + getBound_Start((*(node)->LIS->LB)[0].second) + ");\n";
                    errs() << space + "progress[t_select][" + to_string(currentLoops.size()) << "] = progress[t_select][" + to_string(currentLoops.size()) +  "] % " + getBound_Start((*(node)->LIS->LB)[0].second) + ";\n";
                } else {
                    /* increment the index in the middle level */
                    errs() << space + "progress[t_select][" + to_string(i) + "] = progress[t_select][" + to_string(i) + "] + (progress[t_select][" + to_string(i+1) + "] / " + getBound_Start((*(*(rit-1))->LIS->LB)[0].second) + ");\n";
                    /* module the index in the prev ious level */
                    errs() << space + "progress[t_select][" + to_string(i+1) + "] = progress[t_select][" + to_string(i+1) + "] % " + getBound_Start((*(*(rit-1))->LIS->LB)[0].second) + ";\n";
                }
                i --;
            }
        }
        errs() << "#ifdef DEBUG\n";
        errs() << space + "// cout <<  \"[Thread \" << t_select << \"] next iteration: \";\n";
        errs() << space + "for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {\n";
        errs() << space + "    // cout << *it << \" \";\n";
        errs() << space + "}\n";
        errs() << space + "// cout << endl;\n";
        errs() << "#endif\n";
        // check if is out of bound
        if (allLevelInc) {
            errs() << space + "if (progress[t_select][0] > BLIST[t_select][1]) {\n";
        } else {
            errs() << space + "if (progress[t_select][" << currentLoops.size()<< "]";
            if ((*node->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLE || (*node->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULE) {
                    errs() << " > ";
                } else if ((*node->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGE || (*node->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGE) {
                    errs() << " < ";
                } else if ((*node->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SLT || (*node->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_ULT) {
                    errs() << " >= ";
                } else if ((*node->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_SGT || (*node->LIS->PREDICATE)[0] == llvm::CmpInst::ICMP_UGT) {
                    errs() << " <= ";
                } 
            errs() << getBound_Start((*(node)->LIS->LB)[0].second) << ") {\n";
        }
        errs() << space + "    // remove t_select from the thread pool\n";
        errs() << space + "    candidate_thread_pool_" << node->LoopLevel << ".erase(remove(candidate_thread_pool_" << node->LoopLevel << ".begin(), candidate_thread_pool_" << node->LoopLevel << ".end(), t_select), candidate_thread_pool_" << node->LoopLevel << ".end());\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space + "    cout << \"Remove thread \" << t_select << \" from the candidate thread pool\" << endl;\n";
        errs() << "#endif\n";
        errs() << space + "}\n";
        return;
    }
    
    

    /* Generating Reuse Search */
    bool AllLevelRICodeGen_ref::refRTSearchGen(
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree,
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, 
        string space
                                               ) {
        
        if (!LoopRefTree->isThreadNode) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                hasTNLoop = hasTNLoop || (*it)->isThreadNode;
            }
            if (currentLoops.size() == 0 && LoopRefTree->)
                errs() << "#ifdef DEBUG\n";
                errs() << space + "// cout << \"Count: \" << cnt << endl;\n";
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
                /* Outer-most loops will generate chunk iterations */
                errs() << space + "/* Generating thread local iteration space mapping code */\n";
                errs() << space + "for (int cid = 0; cid < chunk_num; cid++) {\n";
                errs() << space + "    /* Computes bound express for each thread */\n";
                errs() << space + "    for (int t = 0; t < THREAD_NUM; ++t) {\n";
                errs() << space + "        BLIST[t][0] =  " + getBound((*LoopRefTree->LIS->LB)[0].first) + "+ (cid * THREAD_NUM + t) * chunk_size;\n";
                errs() << space + "        BLIST[t][1] = min(" + getBound((*LoopRefTree->LIS->LB)[0].first) + " + (cid * THREAD_NUM + t + 1) * chunk_size, " + getBound((*LoopRefTree->LIS->LB)[0].second) + ") - 1;\n";
                errs() << "#ifdef DEBUG\n";
                errs() << space + "        // cout << \"[Thread \" << t << \"], \" << \"(\" << BLIST[t][0] << \", \"<< BLIST[t][1] << \")\" << endl;\n";
                errs() << "#endif\n";
                errs() << space + "    }\n";
                errs() << space + "    vector<int> thread_pool;\n";
                errs() << space + "    map<int, vector<int>> progress;\n";
            }
            if (hasTNLoop) {
                errs() << space + "/* Generate the Random Interleaving process */\n";
                // if it's not the first loop
                if (find(outloops.begin(), outloops.end(), LoopRefTree) != outloops.end() || isFirstChildLoop) {
                    errs() << space + "vector<int> candidate_thread_pool_" << LoopRefTree->LoopLevel << ";\n";
                }
                errs() << space + "for (int tid = 0; tid < THREAD_NUM; tid++) {\n";
                errs() << space + "    candidate_thread_pool_" << LoopRefTree->LoopLevel << ".push_back(tid);\n";
                errs() << space + "    /* init the progress vector for each thread (" << currentLoops.size() << ") */\n";
                errs() << space + "    progress[tid]";
                string tmp = "";
                if (!containsFirstTN) {
                    tmp += " = {";
                    for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        if (find(outloops.begin(), outloops.end(), *it) != outloops.end()) {
                            tmp += "cid * (THREAD_NUM * chunk_size) + " + getBound((*(*it)->LIS->LB)[0].first) + " + chunk_size * tid";
                        } else {
                            tmp += getBound((*(*it)->LIS->LB)[0].first);
                        }
                        tmp += ", ";
                    }
                    if (currentLoops.size() == 0) {
                        tmp += "cid * (THREAD_NUM * chunk_size) + " + getBound((*LoopRefTree->LIS->LB)[0].first) + " + chunk_size * tid";
                    } else {
                        tmp += getBound((*(LoopRefTree)->LIS->LB)[0].first);
                    }
                    tmp += " };\n";
                } else {
                    tmp += ".push_back(" + getBound((*(LoopRefTree)->LIS->LB)[0].first) + ");\n";
                }
                errs() << tmp;
                errs() << space + "}\n";
                errs() << space + "while ( !candidate_thread_pool_" << LoopRefTree->LoopLevel << ".empty()) {\n";
            }
        } else { // Thread Node
            // find reuses for each ref node
            for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                if ((*it)->AA != NULL) {
                    if (arrayName[(*it)->AA] == refName) {
                        errs() << space + "for(vector<int>::iterator it = candidate_thread_pool_" << LoopRefTree->LoopLevel << ".begin(); it != candidate_thread_pool_" << LoopRefTree->LoopLevel << ".end(); ++it) {\n";
                        errs() << space + "        thread_pool.push_back(*it);\n";
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "    cout << \"[\" << *it << \"] Iteration \" << ";
                        for (unsigned i = 0; i < currentLoops.size(); i++) {
                                errs() << "progress[*it][" + std::to_string(i) + "] << \" \" << ";
                            }
                        errs() << " endl;\n";
                        errs() << "#endif\n";
                        errs() << space + "}\n";
                        errs() << space + "while ( !thread_pool.empty()) {\n";
                        errs() << space + "    int t_select = thread_pool[rand() % thread_pool.size()];\n";
                        errs() << space + "    cnt++;\n";
                        errs() << space + "    access = calAddr" + arrayName[(*it)->AA] + std::to_string(refNumber[(*it)->AA]) + "( ";
                        string tmp = "";
                        for (unsigned i = 0; i < currentLoops.size(); i++) {
                            tmp += "progress[t_select][" + std::to_string(i) + "], ";
                        }
                        if (currentLoops.size() != 0) {
                            tmp.pop_back();
                            tmp.pop_back();
                        }
                        errs() << tmp + ");\n";
                        errs() << space + "    if (LAT.find(access) != LAT.end()) {\n";
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "        cout << \"[REUSE of Addr \" << access << \"] \" << cnt - get<0>(LAT[access]) << \" find @(\" ";
                        for (unsigned i = 0; i < currentLoops.size(); i++) {
                            errs() << "<< progress[t_select][" + std::to_string(i) + "] << \" \"";
                        }
                        errs() << "<< \")\" << endl;\n";
                        errs() << "#endif\n";
                        errs() << space + "        subBlkRT(RT, cnt - get<0>(LAT[access]));\n";
                        errs() << space + "    }\n";
                        errs() << space + "    LAT[access] = make_tuple(cnt, cid);\n";
                         /* remove from thread pool */
                        errs() << space + "    thread_pool.erase(remove(thread_pool.begin(), thread_pool.end(), t_select), thread_pool.end());\n";
                        errs() << "#ifdef DEBUG\n";
                        errs() << space + "    // cout << \"Remove thread \" << t_select << \" from the pool\" << endl;\n";
                        errs() << "#endif\n";
                        errs() << space + "}\n";
                    } else {
                        errs() << space + "cnt += THREAD_NUM;\n";
                    }
                }
            }
        }
        
        
        
        return GenFlag;
    }


    //void AllLevelRICodeGen_ref::refRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, string refName, int useID) {
    void AllLevelRICodeGen_ref::refRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
       // std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops = findLoops(LoopRefTree, refName, useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(), false);
        
        //for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(); lit != loops.end(); lit++) {
        //    errs() << "/* for (" << indvName[(*(*lit)->LIS->IDV)[0]] << ", " << getBound_Start((*(*lit)->LIS->LB)[0].first) << ", " << getBound_Start((*(*lit)->LIS->LB)[0].second)  << ") */\n";
        //}
        
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> sampleIDVs;
        
        string space = "    ";
  
        errs() << space + "uint64_t cnt = 0;\n";
        errs() << "\n";
        errs() << space + "/* Variable used to compute thread-local iteration space (out-most-loops) */\n";
        errs() << space + "auto BLIST = new int[THREAD_NUM][2];\n";
        errs() << space + "/* Vector that contains the interleaved iteration, avoid duplicate declaration */\n";
        errs() << space + "int chunk_size, chunk_num;\n";
        errs() << space + "uint64_t access;\n";
        refRTSearchGen(LoopRefTree, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(), "    ");
        return;
    }
    
    void AllLevelRICodeGen_ref::refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        errs() << "void interleaving() {\n";
        refRTBodyGen(LoopRefTree);
        errs() << "}\n";
        
        return;
    }
    
    void AllLevelRICodeGen_ref::mainGen() {
    
        errs() << "int main() {\n";
        
        string space = "    ";
        errs() << "    interleaving();\n";

#ifdef DumpRTMR
#ifdef REFERENCE_GROUP
        
        errs() << "    refRTDump();\n";
        // errs() << "    gaussian_distr();\n";
        // errs() << "    uniform_distr();\n";
        errs() << "    refRTDump();\n";  
        errs() << "    rtMerge();\n";
#endif
        errs() << "    rtDump();\n";
        errs() << "    RTtoMR_AET();\n";
        errs() << "    dumpMR();\n";
#elif defined(DumpRefLease)
        errs() << "    RL_main(0);\n";
#endif
        
        errs() << "    return 0;\n";
        errs() << "}\n";
        
        return;
    }
    
    
    bool AllLevelRICodeGen_ref::runOnFunction(Function &F) {
    
        errs() << " // Start to generating Static Sampling Code (reference based)\n";
        
        /* reading info from previous passes */
        arrayName = getAnalysis<idxAnalysis::IndexAnalysis>().arrayName;
        arrayExpression = getAnalysis<idxAnalysis::IndexAnalysis>().arrayExpression;
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
        GaussianDistrGen();
        UniformDistrGen();
        rtMergeGen();
#endif
        
#ifdef DumpRTMR
        /* generate rtToMR function */
        rtToMRGen();
        
        /* generate rtDump function */
        rtDumpGen();
        
        /* generate mrDump function */
        mrDumpGen();
#endif
        
        /* generate addr cal function*/
        addrCalFuncGenTop(LoopRefTree);
        
        /* generate rtGen */
        refRTGen(LoopRefTree);
        
        /* generate main function */
        mainGen();

        return false;
    }
    
    void AllLevelRICodeGen_ref::getAnalysisUsage(AnalysisUsage &AU) const {
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
