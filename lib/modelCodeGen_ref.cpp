#include "modelCodeGen_ref.hpp"


using namespace std;

namespace modelCodeGen_ref {

    char ModelCodeGen_ref::ID = 0;
    static RegisterPass<ModelCodeGen_ref> X("modelCodeGen_ref", "static sampling code generating pass (reference based)", false, false);

    ModelCodeGen_ref::ModelCodeGen_ref() : FunctionPass(ID) {}

    void ModelCodeGen_ref::numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
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
    
    void ModelCodeGen_ref::numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
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
    
    void ModelCodeGen_ref::initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
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
            for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                initIndvName(*it);
            }
        }
        
        return;
    }
    
    void ModelCodeGen_ref::initArrayName() {
        
        for ( map<Instruction*,  string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            it->second.replace( find(it->second.begin(), it->second.end(), '.'),  find(it->second.begin(), it->second.end(), '.') +1, 1, '_');
        }
        
        return;
    }
    
    void ModelCodeGen_ref::addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree,  vector<string> indvs) {
        
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
    
    void ModelCodeGen_ref::addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        vector<string> indvs;
        addrCalFuncGen(LoopRefTree, indvs);
        
        return;
    }
    
    void ModelCodeGen_ref::headerGen() {
        
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
        errs() << " map<uint64_t, double> RT;\n";
        errs() << " map<uint64_t, double> MR;\n";
#ifdef REFERENCE_GROUP
        errs() << " map<string, map<uint64_t, double>> refRT;\n";
#endif
        return;
    }

#ifdef DumpRTMR
    void ModelCodeGen_ref::rtHistoGen() {
        
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
#endif

     /* Generate the function to calculate the bins */
    void ModelCodeGen_ref::subBlkRTGen() {
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
    void ModelCodeGen_ref::rtDumpGen() {
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
    
    void ModelCodeGen_ref::rtToMRGen() {
        
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
    
    void ModelCodeGen_ref::mrDumpGen() {
        
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
#endif

#if defined(REFERENCE_GROUP)
    void ModelCodeGen_ref::rtMergeGen() {
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

    string ModelCodeGen_ref::getBound(Value *bound) {
        
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
    
    string ModelCodeGen_ref::getBound_Start(Value *bound) {
        
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

    string ModelCodeGen_ref::getLoopInc(Value *inc) {
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
    
    vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> ModelCodeGen_ref::findLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, string refName, int useID, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops) {
        
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
    void ModelCodeGen_ref::searchReuseDifferentLoopsUpdateFuncGen() {
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
    
    void ModelCodeGen_ref::searchReuseDifferentLoopsCalFuncGen() {
        
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
    
    bool ModelCodeGen_ref::searchReuseDifferentLoopsInitGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, bool GenFlag,  string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> currentLoops, string space) {
        
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
    
    bool ModelCodeGen_ref::searchReuseDifferentLoopsBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag,  string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> currentLoops, string space) {
        
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
    void ModelCodeGen_ref::searchReuseSameLoopInitGen( string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, string space) {
        
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
    
    void ModelCodeGen_ref::searchReuseSameLoopBodyGen( string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops,  string space) {
        
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
    bool ModelCodeGen_ref::refRTSearchGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag, string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, string space) {
        
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
                    if (arrayTypeMap[refName + std::to_string(useID)] != REGULAR) {
                        is_normal_ref = "false";
                    } 
                    errs() << space + "            uint64_t parallel_rt = parallel_predict((" << indvName[(*(*currentLoops.begin())->LIS->IDV)[0]] << "_Start -" << getBound_Start((*(*currentLoops.begin())->LIS->LB)[0].first) << "), " << indvName[(*(*currentLoops.begin())->LIS->IDV)[0]] << ", cnt, " << to_string(outMostLoopPerIterationSpace[refName +  to_string(useID)]) << ", " << to_string(outMostLoopPerIterationSpace[refName +  to_string(refNumber[LoopRefTree->AA])]) << ", "<< is_normal_ref << ");\n";
#ifdef REFERENCE_GROUP
                    // errs() << space + "            refSubBlkRT(refRT, cnt, 1.0, \"" + refName + std::to_string(useID) + "\");\n";
                    errs() << space + "            refRTHistoCal(refRT, parallel_rt, 1.0, \"" + refName + std::to_string(useID) + "\"";
#else
                    // errs() << space + "            subBlkRT(RT, cnt, 1.0);\n";
                    errs() << space + "            rtHistoCal(RT, parallel_rt, 1.0";

#endif
                    // errs() << space + "        subBlkRT(RT, cnt, 1.0";
                    errs() << ");\n";
                    errs() << "#ifdef DEBUG\n";
                    // errs() << space + "            if (parallel_rt != 1 && parallel_rt != 13) {\n";
                    errs() << space + "                cout << \"[\" << parallel_rt << \"] (\" << ";
                    for ( vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            errs() << indvName[(*(*it)->LIS->IDV)[i]] << "_Start";
                            if (it != currentLoops.end()-1) {
                                errs() << "<< \", \" << ";
                            }
                        }
                    }
                    errs() << "<< \") -- (\" << ";
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


    void ModelCodeGen_ref::refRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, string refName, int useID) {
        
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
    
    void ModelCodeGen_ref::refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
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
    
    void ModelCodeGen_ref::mainGen() {
    
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
#ifdef REFERENCE_GROUP
        errs() << "    rtMerge();\n";
        errs() << "    refRTDump();\n";
#endif
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

    void ModelCodeGen_ref::defineRefs(Instruction* arrayInstr) {
        errs() << "/* Array " << arrayName[arrayInstr] << "\t";
        for(vector<string>::iterator vit = arrayAccessVariable[arrayInstr].begin(); vit != arrayAccessVariable[arrayInstr].end(); ++vit) {
            errs() << *vit << " ";
        }
        errs() << "*/ \n";
        for(vector<string>::iterator vit = arrayAccessVariable[arrayInstr].begin(); vit != arrayAccessVariable[arrayInstr].end(); ++vit) {
            // iterate all outermost loops
            for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = outloops.begin(); it != outloops.end(); ++it) {
                // contains outermost loop induction variable
                if (indvName[(*(*it)->LIS->IDV)[0]] == *vit ) {
                    // contains only 1 induction 
                    if (arrayAccessVariable[arrayInstr].size() == 1) {
                        arrayTypeMap[arrayName[arrayInstr] + to_string(refNumber[arrayInstr])] = OUTERMOST_ONLY;
                        return;
                    } else if (vit != arrayAccessVariable[arrayInstr].begin()) {
                        arrayTypeMap[arrayName[arrayInstr] + to_string(refNumber[arrayInstr])] = REVERSE;
                        return;
                    }
                    arrayTypeMap[arrayName[arrayInstr] + to_string(refNumber[arrayInstr])] = REGULAR;
                    return;
                }
            }
        }
        arrayTypeMap[arrayName[arrayInstr] + to_string(refNumber[arrayInstr])] = OUTERMOST_INVARIANT;
        return;
    }
    
    void ModelCodeGen_ref::parallelModelCodeGen() {
        string space = "    ";
        errs() << "int getChunkID(uint64_t i) {\n";
        errs() << space << "return floor(i / (CHUNK_SIZE * THREAD_NUM));\n";
        errs() << "}\n";

        errs() << "int getThreadID(uint64_t i) {\n";
        errs() << space << "return i / CHUNK_SIZE - floor(i / (CHUNK_SIZE * THREAD_NUM))*THREAD_NUM ;\n";
        errs() << "}\n";

        errs() << "int getThreadLocalPos(uint64_t i) {\n";
        errs() << space << "return i \% CHUNK_SIZE;\n";
        errs() << "}\n";

        errs() << "uint64_t parallel_predict(uint64_t i_src, uint64_t i_sink, uint64_t rt, uint64_t lsrc, uint64_t lsink, bool is_normal_ref) {\n";
        errs() << space << "uint64_t parallel_rt = rt;\n";
        errs() << space << "int tsrc = getThreadID(i_src);\n";
        errs() << space << "int tsink = getThreadID(i_sink);\n";
        errs() << space << "int dT = tsink - tsrc;\n";
        errs() << space << "if (!is_normal_ref) {\n";
        errs() << space << space << "if (tsrc < THREAD_NUM - 1) {\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Neighboring Effect\" << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "return 1;\n";
        errs() << space << space << "}\n";
        errs() << space << space << "/*\n";
        errs() << space << space << "} else {\n";
        errs() << space << space << space << "return (rt - 1) * THREAD_NUM + 1;\n";
        errs() << space << space << "}\n";
        errs() << space << space << "*/\n";
        errs() << space << "}\n";
        errs() << space << "/* intra chunk reuse */\n";
        errs() << space << "if (getChunkID(i_src) == getChunkID(i_sink)) {\n"; 
        errs() << space << space << "/* same thread -- scaling effect */\n";
        errs() << space << space <<"if (dT == 0) {\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Scaling Effect\" << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "parallel_rt = rt * THREAD_NUM;\n";
        errs() << space << space << "} else if (getThreadLocalPos(i_src) <= getThreadLocalPos(i_sink)) { // src-sink order\n";
        errs() << space << space << space << "if ((rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + dT) < 0) { printf(\"NORMAL ORDER NEGATIVE PRI\\n\"); }\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Src-Sink Order Folding Effect\" << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "parallel_rt = rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + abs(dT);\n";
        errs() << space << space << "} else { // sink-src order\n";
        errs() << space << space << space << "if ((rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + dT) < 0) { printf(\"REVERSE ORDER NEGATIVE PRI\\n\"); }\n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Sink-Src Order Folding Effect\" << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "parallel_rt = CHUNK_SIZE * lsrc * THREAD_NUM * dT - (rt * THREAD_NUM) - abs(dT);\n";
        errs() << space << space << "}\n";
        errs() << space << "} else { // inter chunk reuse \n";
        errs() << "#ifdef DEBUG\n";
        errs() << space << space << space << "cout << \"Inter Chunk Reuse\" << endl;\n";
        errs() << "#endif\n";
        errs() << space << space << space << "parallel_rt = rt * THREAD_NUM - CHUNK_SIZE * THREAD_NUM * (lsrc*(THREAD_NUM - tsrc) + lsink * tsink) + CHUNK_SIZE * THREAD_NUM * lsink + dT;\n";
        errs() << space << "}\n";
        errs() << space << "return parallel_rt;\n";
        errs() << "}\n";
    }
    
    bool ModelCodeGen_ref::runOnFunction(Function &F) {
    
        errs() << " // Start to generating Static Sampling Code (reference based)\n";
        
        /* reading info from previous passes */
        arrayName = getAnalysis<idxAnalysis::IndexAnalysis>().arrayName;
        arrayExpression = getAnalysis<idxAnalysis::IndexAnalysis>().arrayExpression;
        arrayAccessVariable= getAnalysis<idxAnalysis::IndexAnalysis>().arrayAccessVariable;
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree;
        sampleNum = getAnalysis<sampleNumAnalysis::SampleNumberAnalysis>().sampleNum;
#ifdef PARALLEL
        map<Instruction*, uint64_t> tempL = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().outMostLoopPerIterationSpace;;
        outloops = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().outMostLoops;
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

        parallelModelCodeGen();

        /* generate rtHistoCal function */
        rtHistoGen();

        /* generate subBlkRT function */
        subBlkRTGen();

#if defined(REFERENCE_GROUP)
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
    
    void ModelCodeGen_ref::getAnalysisUsage(AnalysisUsage &AU) const {
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
