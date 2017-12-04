#include "ssCodeGen_ref.hpp"

namespace ssCodeGen_ref {

    char StaticSamplingCodeGen_ref::ID = 0;
    static RegisterPass<StaticSamplingCodeGen_ref> X("ssCodeGen_ref", "static sampling code generating pass (reference based)", false, false);

    StaticSamplingCodeGen_ref::StaticSamplingCodeGen_ref() : FunctionPass(ID) {}

    void StaticSamplingCodeGen_ref::numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
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
    
    void StaticSamplingCodeGen_ref::numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
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
    
    void StaticSamplingCodeGen_ref::initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree == NULL) {
            return;
        }
        
        if (LoopRefTree->L != NULL) {
            for (std::vector<Value *>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvName[*it] = (*it)->getName();
                indvName[*it].erase(std::find(indvName[*it].begin(), indvName[*it].end(), '.'));
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                initIndvName(*it);
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::initArrayName() {
        
        for (std::map<Instruction*, std::string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            it->second.replace(std::find(it->second.begin(), it->second.end(), '.'), std::find(it->second.begin(), it->second.end(), '.') +1, 1, '_');
        }
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<string> indvs) {
        
        if (LoopRefTree->L != NULL) {
            for (std::vector<Value*>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvs.push_back(indvName[(*it)]);
            }
        }
        
        if (LoopRefTree->AA != NULL) {
            errs() << "int calAddr" + arrayName[LoopRefTree->AA] + std::to_string(refNumber[LoopRefTree->AA]) + "( ";
            string arguments = "";
            for (std::vector<string>::iterator it = indvs.begin(), eit = indvs.end(); it != eit; ++it) {
                arguments += "int " + (*it) + ", ";
            }
            arguments.pop_back();
            arguments.pop_back();
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
    
    void StaticSamplingCodeGen_ref::addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        std::vector<string> indvs;
        addrCalFuncGen(LoopRefTree, indvs);
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::headerGen() {
        
        errs() << "#include <map>\n";
        errs() << "#include <set>\n";
        errs() << "#include <cstdlib>\n";
        errs() << "#include <iostream>\n";
#ifdef PARALLEL_CXX_THREAD
        errs() << "#include <thread>\n";
        errs() << "#include <mutex>\n";
#endif
        errs() << "using namespace std;\n";
        
#ifdef PARALLEL_CXX_THREAD
        errs() << "std::mutex mtx;\n";
#endif
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::rtHistoGen() {
        
        errs() << "std::map<uint64_t, double> RT;\n";
        errs() << "std::map<uint64_t, double> MR;\n";
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
        errs() << "    std::unique_lock<std::mutex> lck (mtx,std::defer_lock);\n";
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
        
        return;
    }
    
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
        errs() << space + "for (uint64_t c = 0; c <= max_RT; c++) {\n";
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
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(bound)) {
            return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
        }
        return "";
    }
    
    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> StaticSamplingCodeGen_ref::findLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops) {
        
        if (LoopRefTree->L != NULL) {
            loops.push_back(LoopRefTree);
        }
        if (LoopRefTree->AA != NULL) {
            if (arrayName[LoopRefTree->AA] == refName && refNumber[LoopRefTree->AA] == useID) {
                return loops;
            } else {
                return std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>();
            }
        }
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopRes;
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loopTmp = findLoops(*it, refName, useID, loops);
                if (loopTmp.size() != 0) {
                    loopRes = loopTmp;
                }
            }
        }
        return loopRes;
    }
    
    
    /* Search result reuse */
    void StaticSamplingCodeGen_ref::searchReuseSameLoopInitGen(std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, string space) {
        
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator ait = loops.back()->next->begin(), eit = loops.back()->next->end(); ait != eit; ++ait) {
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
    
    void StaticSamplingCodeGen_ref::searchReuseSameLoopBodyGen(std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops,  string space) {
        
        
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator ait = loops.back()->next->begin(), eit = loops.back()->next->end(); ait != eit; ++ait) {
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
                    
                    errs() << space + "        rtHistoCal(prev_cnt_" + refName + std::to_string(refNumber[(*ait)->AA]) + ");\n";
//                    errs() << space + "        cout << \"Find match\\n\";\n";
                    errs() << space + "        goto EndSample;\n";
                    
                    errs() << space + "    }\n";
                    errs() << space + "}\n";
                }
            }
        }
        
        return;
    }
    
    
    /* Generating Reuse Search */
    bool StaticSamplingCodeGen_ref::refRTSearchGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag, string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, string space) {
        
        if (loops.size() != 0 && GenFlag == false) {
            if (LoopRefTree == loops[0]) {
                GenFlag = true;
            }
        }
        
        if (GenFlag == true) {
            
            /* Generate loop */
            if (LoopRefTree->L != NULL) {
                
                errs() << space + "{\n";
                
                string loopNum = std::to_string(loopNumber[LoopRefTree->L]);
                
                if (std::find(loops.begin(), loops.end(), LoopRefTree) != loops.end()) {
                    if (LoopRefTree == loops[0]) {
                        errs() << space + "int " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = ";
                        errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] + "_Start;\n";
                    } else {
                        errs() << space + "int " + indvName[(*LoopRefTree->LIS->IDV)[0]] + "LB" + loopNum + " = ";
                        errs() << getBound((*LoopRefTree->LIS->LB)[0].first) + ";\n";
                        errs() << space + "if ( ";
                        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
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
                errs() << " < ";
                errs() << getBound((*LoopRefTree->LIS->LB)[0].second);
                errs() << "; ";
                errs() << indvName[(*LoopRefTree->LIS->IDV)[0]] + "++";
                errs() << ") {\n";
            }
            
            /* Generate addr checking instruction */
            if (LoopRefTree->AA != NULL) {
                
                if (arrayName[LoopRefTree->AA] == refName) {
                
                    
                    errs() << space + "if (cntStart == true) {\n";
                    errs() << space + "    cnt++;\n";
                    
                    errs() << space + "    if ( calAddr" + refName + std::to_string(refNumber[LoopRefTree->AA]) + "( ";
                    string tmp = "";
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = currentLoops.begin(), eit = currentLoops.end(); it != eit; ++it) {
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
                    
                    errs() << space + "        rtHistoCal(cnt);\n";
#ifdef SEARCH_REUSE
                    if (std::find(loops.back()->next->begin(), loops.back()->next->end(), LoopRefTree) != loops.back()->next->end()) {
                        errs() << space + "prev_cnt_" + refName + std::to_string(refNumber[LoopRefTree->AA]) + " = cnt;\n";
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
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
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
    
    void StaticSamplingCodeGen_ref::refRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, string refName, int useID) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops = findLoops(LoopRefTree, refName, useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>());
        
#ifdef SEARCH_REUSE
        errs() << "    /* Generating search reuse init code */\n";
        searchReuseSameLoopInitGen(refName, useID, loops, "    ");
#endif
        
#if SAMPLING ==2
        errs() << "    /* Generating sampling loop */\n";
        string space = "    ";
        errs() << space + "set<string> record;\n";
        errs() << space + "for ( int s = 0; s < ";
        int sampling_num = 1;
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                sampling_num *= (int) ((stoi(getBound((*(*lit)->LIS->LB)[i].second)) - stoi(getBound((*(*lit)->LIS->LB)[i].first))) * RANDOM_REF_SAMPLING_RATE);
            }
        }
        errs() << std::to_string(sampling_num);
        errs() << ";) {\n";
        
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                errs() << space + "    int " + indvName[(*(*lit)->LIS->IDV)[i]] + "_Start" + " = ";
                errs() << "rand() % (" + (getBound((*(*lit)->LIS->LB)[i].second) + " - " + getBound((*(*lit)->LIS->LB)[i].first)) + ") + " + getBound((*(*lit)->LIS->LB)[i].first);
                errs() << ";\n";
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
        errs() << space + "while ( record.find(idx_string) != record.end() ) {\n";
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
            for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                errs() << space + "    " + indvName[(*(*lit)->LIS->IDV)[i]] + "_Start" + " = ";
                errs() << "rand() % (" + (getBound((*(*lit)->LIS->LB)[i].second) + " - " + getBound((*(*lit)->LIS->LB)[i].first)) + ") + " + getBound((*(*lit)->LIS->LB)[i].first);
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
        errs() << space + "record.insert( idx_string );\n";
#endif
  
        errs() << space + "uint64_t cnt = 0;\n";
        errs() << space + "bool cntStart = false;\n";
        errs() << "\n";
        
#ifdef SEARCH_REUSE
        errs() << "        /* Generating search reuse body code */\n";
        searchReuseSameLoopBodyGen(refName, useID, loops, "        ");
#endif
        
        errs() << space + "/* Generating reuse search code */\n";
        errs() << "\n";
        refRTSearchGen(LoopRefTree, false, refName, useID, loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>(), "    ");
        
        errs() << "EndSample:\n";
        errs() << space + "s++;\n";
        errs() << space + "}\n";
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        for (std::map<string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < it->second; i++) {
                errs() << "void ref_" + it->first + std::to_string(i) + "() {\n";
                refRTBodyGen(LoopRefTree, it->first, i);
                errs() << "}\n";
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen_ref::mainGen() {
    
        errs() << "int main() {\n";
        
        string space = "    ";
        
        for (std::map<string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < (*it).second; i++) {
                errs() << space + "ref_" + (*it).first + std::to_string(i) + "();\n";
            }
        }
        
        errs() << "    rtDump();\n";
        errs() << "    RTtoMR_AET();\n";
        errs() << "    dumpMR();\n";
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
        
        /* init */
        initArrayName();
        numberRefToSameArray(LoopRefTree);
        numberLoops(LoopRefTree);
        initIndvName(LoopRefTree);
        
        /* generate headers */
        headerGen();
        
        /* generate rtHistoCal function */
        rtHistoGen();
        
        /* generate rtToMR function */
        rtToMRGen();
        
        /* generate rtDump function */
        rtDumpGen();
        
        /* generate mrDump function */
        mrDumpGen();
        
        /* generate addr cal function*/
        addrCalFuncGenTop(LoopRefTree);
        
        /* generate rtGen */
        refRTGen(LoopRefTree);
        
        /* generate main function */
        mainGen();
        
        return false;
    }
    
    void StaticSamplingCodeGen_ref::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
        
        return;
    }
    
}
