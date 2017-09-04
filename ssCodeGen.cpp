//
//  ssCodeGen.cpp
//  LLVM
//
//  Created by dongchen on 8/16/17.
//
//

#include "ssCodeGen.hpp"

namespace ssCodeGen {
    
    char StaticSamplingCodeGen::ID = 0;
    static RegisterPass<StaticSamplingCodeGen> X("ssCodeGen", "static sampling code generating pass", false, false);
    
    StaticSamplingCodeGen::StaticSamplingCodeGen() : FunctionPass(ID) {}
    
    string StaticSamplingCodeGen::getBound(Value *bound) {
        
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
    
    void StaticSamplingCodeGen::useLoopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops) {
        
        if (LoopRefTree->L != NULL) {
                loops.push_back(LoopRefTree);
        }
        
        if (LoopRefTree->AA != NULL) {
            if (arrayName[LoopRefTree->AA] == refName && refNumber[LoopRefTree->AA] == useID) {
#if SAMPLING == 0
                std::string space = "    ";
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
                    for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                        errs() << space + "for ( int " + indvName[(*(*lit)->LIS->IDV)[i]] << " = ";
                        errs() << getBound((*(*lit)->LIS->LB)[i].first) << "; ";
                        errs() << indvName[(*(*lit)->LIS->IDV)[i]] << " < " + getBound((*(*lit)->LIS->LB)[i].second) + "; ";
                        errs() << indvName[(*(*lit)->LIS->IDV)[i]] << "++) {\n";
                    }
                    space += "    ";
                }
#elif SAMPLING == 1
                std::string space = "    ";
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
                    for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                        errs() << space + "for ( int " + indvName[(*(*lit)->LIS->IDV)[i]] << " = ";
                        errs() << getBound((*(*lit)->LIS->LB)[i].first) << "; ";
                        errs() << indvName[(*(*lit)->LIS->IDV)[i]] << " < " + getBound((*(*lit)->LIS->LB)[i].second) + "; ";
                        errs() << indvName[(*(*lit)->LIS->IDV)[i]];
                        errs() << " += ";
                        errs() << "1 /" +  std::to_string(UNIFORM_SAMPLING_RATE);
                        errs() << ") {\n";
                    }
                    space += "    ";
                }
#elif SAMPLING == 2
                std::string space = "    ";
                errs() << space + "set<string> record;\n";
                errs() << space + "for ( int s = 0; s < ";
                int sampling_num = 1;
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
                    for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                        sampling_num *= (int) ((stoi(getBound((*(*lit)->LIS->LB)[i].second)) - stoi(getBound((*(*lit)->LIS->LB)[i].first))) * RANDOM_SAMPLING_RATE);
                    }
                }
                errs() << std::to_string(sampling_num);
                errs() << "; s++) {\n";
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
                    for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                        errs() << space + "    int " + indvName[(*(*lit)->LIS->IDV)[i]] + " = ";
                        errs() << "rand() % (" + (getBound((*(*lit)->LIS->LB)[i].second) + " - " + getBound((*(*lit)->LIS->LB)[i].first)) + ") + " + getBound((*(*lit)->LIS->LB)[i].first);
                        errs() << ";\n";
                    }
                }
                
                std::string idx_string_tmp = "";
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
                    for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                        idx_string_tmp += "std::to_string(";
                        idx_string_tmp += indvName[(*(*lit)->LIS->IDV)[i]];
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
                        errs() << space + "    " + indvName[(*(*lit)->LIS->IDV)[i]] + " = ";
                        errs() << "rand() % (" + (getBound((*(*lit)->LIS->LB)[i].second) + " - " + getBound((*(*lit)->LIS->LB)[i].first)) + ") + " + getBound((*(*lit)->LIS->LB)[i].first);
                        errs() << ";\n";
                    }
                }
                
                idx_string_tmp = "";
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
                    for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                        idx_string_tmp += "std::to_string(";
                        idx_string_tmp += indvName[(*(*lit)->LIS->IDV)[i]];
                        idx_string_tmp += ") + \"_\" + ";
                    }
                }
                idx_string_tmp.pop_back();
                idx_string_tmp.pop_back();
                
                errs() << space + "    idx_string = " + idx_string_tmp + ";\n";
                errs() << space + "}\n";
                errs() << space + "record.insert( idx_string );\n";
#endif
            }
            return;
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                useLoopGen(*it, refName, useID, loops);
            }
        }
        
        return;
    }
    
    
    loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* StaticSamplingCodeGen::findRef(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int refID) {
        
        if (LoopRefTree == NULL) {
            return NULL;
        }
        
        if (LoopRefTree->AA != NULL) {
            
            if (arrayName[LoopRefTree->AA] == refName && refNumber[LoopRefTree->AA] == refID) {
                return LoopRefTree;
            }
        }
        
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* result = NULL;
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* tmp = findRef(*it, refName, refID);
                if (tmp != NULL) {
                    result = tmp;
                }
            }
        }
        
        return result;
    }
    
    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> StaticSamplingCodeGen::findRefLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int refID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops) {
        
        if (LoopRefTree->L != NULL) {
            loops.push_back(LoopRefTree);
        }
        
        if (LoopRefTree->AA != NULL) {
            
            if (arrayName[LoopRefTree->AA] == refName && refNumber[LoopRefTree->AA] == refID) {
                return loops;
            }
        }
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> result;
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> tmp = findRefLoops(*it, refName, refID, loops);
                if (!tmp.empty()) {
                    result = tmp;
                }
            }
        }
        
        return result;
    }
    
    bool StaticSamplingCodeGen::reuseLoopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int useID, int reuseID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops = findRefLoops(LoopRefTree, refName, useID, loops);
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops = findRefLoops(LoopRefTree, refName, reuseID, loops);
        
        bool breakflag = false;
        if (useLoops.size() != 0 && reuseLoops.size() != 0 && loopOrder[useLoops[0]] > loopOrder[reuseLoops[0]]) {
            breakflag = true;
        }
        
        std::string space = "    ";
#if (SAMPLING == 0 || SAMPLING == 1)
        for (unsigned long i = 0; i < useLoops.size(); i++) {
            space += "    ";
        }
#elif SAMPLING == 2
        space += "    ";
#endif

        if (breakflag == true) {
#if (SAMPLING == 0 || SAMPLING == 1)
            for (unsigned long i = 0; i < useLoops.size(); i++) {
                errs() << "    break;\n";
                errs() << "}\n";
                space.pop_back();
                space.pop_back();
                space.pop_back();
                space.pop_back();
            }
#elif SAMPLING == 2
            errs() << "    break;\n";
            errs() << "}\n";
#endif
            return false;
        }
        
        errs() << space + "bool findReuseFlag = false;\n";
        
        std::map<Value*, std::string> loopLB;
        std::map<Value*, std::string> loopUB;
        
        for (unsigned long i = 0; i < reuseLoops.size(); i++) {
            
            if (useLoops[i] == reuseLoops[i]) {
                
                if (i == 0) {
                    for (unsigned long j = 0; j < reuseLoops[i]->LIS->IDV->size(); j++) {
                        errs() << space + "for ( int " + indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                        errs() << " = ";
                        
                        std::string tmpLB = indvName[(*useLoops[i]->LIS->IDV)[j]];

                        if (i + 1 != reuseLoops.size()) {
                            if (useLoops[i+1] != reuseLoops[i+1] && refTotalOrder[findRef(LoopRefTree, refName, useID)] >= refTotalOrder[findRef(LoopRefTree, refName, reuseID)]) {
                                tmpLB += " + 1";
                            }
                        } else if (i + 1 == reuseLoops.size()) {
                            if (refTotalOrder[findRef(LoopRefTree, refName, useID)] >= refTotalOrder[findRef(LoopRefTree, refName, reuseID)]) {
                                tmpLB += " + 1";
                            }
                        }
                        errs() << tmpLB;
                        loopLB[(*reuseLoops[i]->LIS->IDV)[j]] = tmpLB;
                        
                        errs() << "; ";
                        errs() << indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                        errs() << " < ";
                        errs() << getBound((*reuseLoops[i]->LIS->LB)[j].second);
                     
                        errs() << "; ";
                        errs() << indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse" + "++";
                        errs() << ") {\n";
                        space += "    ";
                    }
                }
                
                if (i > 0) {
                    
                    for (unsigned long j = 0; j < reuseLoops[i]->LIS->IDV->size(); j++) {
                        errs() << space + "int " + indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                        errs() << ";\n";
                    }
                    
                    errs() << space + "if (";
                    
                    for (unsigned long ii = 0; ii < i; ii++) {
                        for (unsigned long j = 0; j < reuseLoops[ii]->LIS->IDV->size(); j++) {
                            errs() << indvName[(*reuseLoops[ii]->LIS->IDV)[j]] + "reuse" + " == " + loopLB[(*reuseLoops[ii]->LIS->IDV)[j]];
                            if (!(ii == i - 1 && j == reuseLoops[ii]->LIS->IDV->size() - 1)) {
                                errs() << " && ";
                            }
                        }
                    }
                    
                    errs() << ") {\n";
                    
                    for (unsigned long j = 0; j < reuseLoops[i]->LIS->IDV->size(); j++) {
                        errs() << space + "    " + indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                        errs() << " = ";
                        
                        std::string tmpLB = indvName[(*useLoops[i]->LIS->IDV)[j]];
                        
                        if (i + 1 != reuseLoops.size()) {
                            if (useLoops[i+1] != reuseLoops[i+1] && refTotalOrder[findRef(LoopRefTree, refName, useID)] >= refTotalOrder[findRef(LoopRefTree, refName, reuseID)]) {
                                tmpLB += " + 1";
                            }
                        } else if (i + 1 == reuseLoops.size()) {
                            if (refTotalOrder[findRef(LoopRefTree, refName, useID)] >= refTotalOrder[findRef(LoopRefTree, refName, reuseID)]) {
                                tmpLB += " + 1";
                            }
                        }
                        errs() << tmpLB;
                        loopLB[(*useLoops[i]->LIS->IDV)[j]] = tmpLB;
                        
                        errs() << ";\n";
                        
                    }
                    
                    errs() << space + "} else {\n";
                    
                    for (unsigned long j = 0; j < reuseLoops[i]->LIS->IDV->size(); j++) {
                        errs() << space + "    " + indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                        errs() << " = ";
                        errs() << getBound((*reuseLoops[i]->LIS->LB)[j].first);
                        errs() << ";\n";
                    }
                    
                    errs() << space + "}\n";
                    
                    for (unsigned long j = 0; j < reuseLoops[i]->LIS->IDV->size(); j++) {
                        
                        errs() << space + "for ( ";
                        
                        errs() << "; ";
                        errs() << indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                        errs() << " < ";
                        errs() << getBound((*reuseLoops[i]->LIS->LB)[j].second);
                        
                        errs() << "; ";
                        errs() << indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse" + "++";
                        errs() << ") {\n";
                        space += "    ";
                    }
                }
                
            } else {
                for (unsigned long j = 0; j < reuseLoops[i]->LIS->IDV->size(); j++) {
                    errs() << space + "for ( int " + indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                    errs() << " = ";
                    errs() << getBound((*reuseLoops[i]->LIS->LB)[j].first) << "; ";
                    errs() << indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                    errs() << " < ";
                    errs() << getBound((*reuseLoops[i]->LIS->LB)[j].second);
                    errs() << "; ";
                    errs() << indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse" + "++";
                    errs() << ") {\n";
                    space += "    ";
                }
            }
            
        }
        
        return true;
    }
    
    
    void StaticSamplingCodeGen::checkIntervenRefGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops, int useID, int reuseID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTreeRoot) {
        
        if (LoopRefTree == NULL) {
            return;
        }
        
        if (LoopRefTree->L != NULL) {
            loops.push_back(LoopRefTree);
        }
        
        if (LoopRefTree->AA != NULL) {
            if (arrayName[LoopRefTree->AA] == refName && refNumber[LoopRefTree->AA] != reuseID) {
                
                /* need to fix when loop size is not zero */
                /* generate loops */
                std::string space = "    ";
                if (loops.size() != 0 && useLoops.size() != 0 && reuseLoops.size() != 0) {
                    if (loopOrder[loops[0]] >= loopOrder[useLoops[0]] && loopOrder[loops[0]] <= loopOrder[reuseLoops[0]]) {
                        
                        std::map<Value*, std::string> loopLB;
                        std::map<Value*, std::string> loopUB;
                        
                        for (unsigned long i = 0; i < loops.size(); i++) {
                            
                            if (i > 0) {
                                
                                for (unsigned long j = 0; j < loops[i]->LIS->IDV->size(); j++) {
                                    errs() << space + "int ";
                                    errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "Interven;\n";
                                }
                                
                                for (unsigned long j = 0; j < loops[i]->LIS->IDV->size(); j++) {
                                    errs() << space + "int ";
                                    errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "IntervenUB;\n";
                                }
                                
                                
                                /* Generate LB condition */
                                
                                errs() << space + "if (" ;
                            
                                for (unsigned long ii = 0; ii < i; ii++) {
                                    for (unsigned long j = 0; j < loops[ii]->LIS->IDV->size(); j++) {
                                        errs() << indvName[(*loops[ii]->LIS->IDV)[j]] + "Interven" + " == " + loopLB[(*loops[ii]->LIS->IDV)[j]];
                                        if (!(ii == i - 1 && j == loops[ii]->LIS->IDV->size() - 1)) {
                                            errs() << " && ";
                                        }
                                    }
                                }
                                errs() << ") {\n";
                                
                                for (unsigned long j = 0; j < loops[i]->LIS->IDV->size(); j++) {
                                    errs() << space + "    ";
                                    errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "Interven" + " = ";
                                    
                                    std::string intervenLB = "";
                                    
                                    if (i < useLoops.size()) {
                                        if (useLoops[i] == loops[i]) {
                                            intervenLB += indvName[(*useLoops[i]->LIS->IDV)[j]];
                                            
                                            if (i + 1 < useLoops.size() && i + 1 < loops.size()) {
                                                if (useLoops[i+1] != loops[i+1] && refTotalOrder[findRef(LoopRefTreeRoot, refName, useID)] >= refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenLB += " + 1 ";
                                                }
                                            } else if (i + 1 == useLoops.size()) {
                                                if (refTotalOrder[findRef(LoopRefTreeRoot, refName, useID)] >= refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenLB += " + 1 ";
                                                }
                                            }
                                            /* new add 2 */
                                            else if (i + 1 == loops.size()) {
                                                if (refTotalOrder[findRef(LoopRefTreeRoot, refName, useID)] >= refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenLB += " + 1 ";
                                                }
                                            }
                                            
                                        } else {
                                            
                                            intervenLB += getBound((*loops[i]->LIS->LB)[j].first);
                                            
                                        }
                                        
                                    } else {
                                        
                                        intervenLB += getBound((*loops[i]->LIS->LB)[j].first);
                                        
                                    }
                                    
                                    loopLB[(*loops[i]->LIS->IDV)[j]] = intervenLB;
                                    
                                    errs() << intervenLB;
                                    errs() << ";\n";
                                }
                                
                                errs() << space + "} else {\n";
                                for (unsigned long j = 0; j < loops[i]->LIS->IDV->size(); j++) {
                                    errs() << space + "    ";
                                    errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "Interven" + " = ";
                                    errs() << getBound((*loops[i]->LIS->LB)[j].first);
                                    errs() << ";\n";
                                    errs() << space + "}\n";
                                }
                                
                                
                                /* Generate UB condition */
                                
                                errs() << space + "if (" ;
                                
                                for (unsigned long ii = 0; ii < i; ii++) {
                                    for (unsigned long j = 0; j < loops[ii]->LIS->IDV->size(); j++) {
                                        errs() << indvName[(*loops[ii]->LIS->IDV)[j]] + "Interven" + " == " + loopUB[(*loops[ii]->LIS->IDV)[j]];
                                        if (!(ii == i - 1 && j == loops[ii]->LIS->IDV->size() - 1)) {
                                            errs() << " && ";
                                        }
                                    }
                                }
                                errs() << ") {\n";
                                
                                for (unsigned long j = 0; j < loops[i]->LIS->IDV->size(); j++) {
                                    errs() << space + "    ";
                                    errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "IntervenUB" + " = ";
                                    
                                    std::string intervenUB = "";
                                    
                                    if (i < reuseLoops.size()) {
                                        if (reuseLoops[i] == loops[i]) {
                                            intervenUB += indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                                            
                                            if (i + 1 < reuseLoops.size() && i + 1 < loops.size()) {
                                                if (reuseLoops[i+1] != loops[i+1] && refTotalOrder[findRef(LoopRefTreeRoot, refName, reuseID)] < refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenUB += " - 1 ";
                                                }
                                            } else if (i + 1 == reuseLoops.size()) {
                                                
                                                if (refTotalOrder[findRef(LoopRefTreeRoot, refName, reuseID)] < refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenUB += " - 1 ";
                                                }
                                            }
                                            /* new add 2 */
                                            else if (i + 1 == loops.size()) {
                                                
                                                if (refTotalOrder[findRef(LoopRefTreeRoot, refName, reuseID)] < refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenUB += " - 1 ";
                                                }
                                            }
                                            
                                        } else {
                                            
                                            intervenUB += getBound((*loops[i]->LIS->LB)[j].second);
                                            
                                        }
                                    } else {
                                        
                                        intervenUB += getBound((*loops[i]->LIS->LB)[j].second);
                                        
                                    }
                                    loopUB[(*loops[i]->LIS->IDV)[j]] = intervenUB;
                                    errs() << intervenUB;
                                    errs() << ";\n";
                                }
                                
                                errs() << space + "} else {\n";
                                
                                for (unsigned long j = 0; j < loops[i]->LIS->IDV->size(); j++) {
                                    errs() << space + "    ";
                                    errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "IntervenUB" + " = ";
                                    errs() << getBound((*loops[i]->LIS->LB)[j].second) + "- 1";
                                    errs() << ";\n";
                                }
                                
                                errs() << space + "}\n";
                                
                            }
                            
                            errs() << space + "for(";
                            for (unsigned long j = 0; j < loops[i]->LIS->IDV->size(); j++) {
                                
                                if (i == 0) {
                                    errs() << "int ";
                                
                                    errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "Interven" + " = ";
                                
                                    std::string intervenLB = "";
                                
                                    if (i < useLoops.size()) {
                                        if (useLoops[i] == loops[i]) {
                                            intervenLB += indvName[(*useLoops[i]->LIS->IDV)[j]];
                                            
                                            if (i + 1 < useLoops.size() && i + 1 < loops.size()) {
                                                if (useLoops[i+1] != loops[i+1] && refTotalOrder[findRef(LoopRefTreeRoot, refName, useID)] >= refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenLB += " + 1 ";
                                                }
                                            } else if (i + 1 == useLoops.size()) {
                                                if (refTotalOrder[findRef(LoopRefTreeRoot, refName, useID)] >= refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenLB += " + 1 ";
                                                }
                                            }
                                            /* new add */
                                            else if (i + 1 == loops.size()) {
                                                if (refTotalOrder[findRef(LoopRefTreeRoot, refName, useID)] >= refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenLB += " + 1 ";
                                                }
                                            }
                    
                                        } else {
                                            
                                            intervenLB += getBound((*loops[i]->LIS->LB)[j].first);
                                        
                                        }
                                    
                                    } else {
                                    
                                        intervenLB += getBound((*loops[i]->LIS->LB)[j].first);
                                    
                                    }
                                
                                    loopLB[(*loops[i]->LIS->IDV)[j]] = intervenLB;
                                    errs() << intervenLB;
                                }
                                    
                                errs() << "; ";
                                errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "Interven";
                                errs() << " <= ";
                                
                                if (i == 0) {
                                    std::string intervenUB = "";
                                
                                    if (i < reuseLoops.size()) {
                                        if (reuseLoops[i] == loops[i]) {
                                            intervenUB += indvName[(*reuseLoops[i]->LIS->IDV)[j]] + "reuse";
                                        
                                            if (i + 1 < reuseLoops.size() && i + 1 < loops.size()) {
                                                if (reuseLoops[i+1] != loops[i+1] && refTotalOrder[findRef(LoopRefTreeRoot, refName, reuseID)] < refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenUB += " - 1 ";
                                                }
                                            } else if (i + 1 == reuseLoops.size()) {
                                                if (refTotalOrder[findRef(LoopRefTreeRoot, refName, reuseID)] < refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenUB += " - 1 ";
                                                }
                                            }
                                            /* new add 2*/
                                            else if (i + 1 == loops.size()) {
                                                if (refTotalOrder[findRef(LoopRefTreeRoot, refName, reuseID)] < refTotalOrder[findRef(LoopRefTreeRoot, refName, refNumber[LoopRefTree->AA])]) {
                                                    intervenUB += " - 1 ";
                                                }
                                            }
                                        
                                        } else {
                                        
                                            intervenUB += getBound((*loops[i]->LIS->LB)[j].second);
                                        
                                        }
                                    } else {
                                    
                                        intervenUB += getBound((*loops[i]->LIS->LB)[j].second);
                                    
                                    }
                                    loopUB[(*loops[i]->LIS->IDV)[j]] = intervenUB;
                                    errs() << intervenUB;
                                } else {
                                    errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "IntervenUB";
                                }
                                
                                errs() << "; ";
                                errs() << indvName[(*loops[i]->LIS->IDV)[j]] + "Interven" + "++";
                                errs() << ") {\n";
                                space += "    ";
                            }
                        }
                        
                        /* generate check same loc function */
                        
                        errs() << space + "if( calAddr" + refName + std::to_string(refNumber[LoopRefTree->AA]);
                        errs() << "(";
                        std::string tmp = "";
                        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                            for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                                tmp += indvName[(*(*it)->LIS->IDV)[i]];
                                tmp += "Interven";
                                tmp += ", ";
                            }
                        }
                        tmp.pop_back();
                        tmp.pop_back();
                        
                        errs() << tmp + ")";
                        errs() << " == ";
                        
                        errs() << "calAddr" + refName + std::to_string(useID);
                        errs() << "(";
                        tmp = "";
                        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = useLoops.begin(), eit = useLoops.end(); it != eit; ++it) {
                            for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                                tmp += indvName[(*(*it)->LIS->IDV)[i]];
                                tmp += ", ";
                            }
                        }
                        tmp.pop_back();
                        tmp.pop_back();
                        
                        errs() << tmp + ")";
                        errs() << ") {\n";
                        
                        errs() << space + "    return true;\n";
                        errs() << space + "}\n";
                        
                        /* finalize */
                        
                        for (unsigned long i = 0; i < loops.size(); i++){
                            space.pop_back();
                            space.pop_back();
                            space.pop_back();
                            space.pop_back();
                            errs() << space + "}\n";
                        }
                    }
                }
                
                /* generate check same loc function */
                /*
                errs() << space + "if( calAddr" + refName + std::to_string(refNumber[LoopRefTree->AA]);
                errs() << "(";
                std::string tmp = "";
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = loops.begin(), eit = loops.end(); it != eit; ++it) {
                    for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                        tmp += indvName[(*(*it)->LIS->IDV)[i]];
                        tmp += "Interven";
                        tmp += ", ";
                    }
                }
                tmp.pop_back();
                tmp.pop_back();
                
                errs() << tmp + ")";
                errs() << " == ";
                
                errs() << "calAddr" + refName + std::to_string(useID);
                errs() << "(";
                tmp = "";
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = useLoops.begin(), eit = useLoops.end(); it != eit; ++it) {
                    for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                        tmp += indvName[(*(*it)->LIS->IDV)[i]];
                        tmp += ", ";
                    }
                }
                tmp.pop_back();
                tmp.pop_back();
                
                errs() << tmp + ")";
                errs() << ") {\n";
                
                errs() << space + "    return true;\n";
                errs() << space + "}\n";
                */
                /* finalize */
                /*
                for (unsigned long i = 0; i < loops.size(); i++){
                    space.pop_back();
                    space.pop_back();
                    space.pop_back();
                    space.pop_back();
                    errs() << space + "}\n";
                }
                 */
            }
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                checkIntervenRefGen((*it), refName, useLoops, reuseLoops, useID, reuseID, loops, LoopRefTreeRoot);
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::checkIntervenBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int useID, int reuseID) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops;
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops = findRefLoops(LoopRefTree, refName, useID, loops);
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops = findRefLoops(LoopRefTree, refName, reuseID, loops);
        
        checkIntervenRefGen(LoopRefTree, refName, useLoops, reuseLoops, useID, reuseID, loops, LoopRefTree);
        
        return;
    }
    
    void StaticSamplingCodeGen::checkIntervenTopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
    
        for (std::map<std::string, int>::iterator ref_it = refToSameArrayCnt.begin(), ref_eit = refToSameArrayCnt.end(); ref_it != ref_eit; ++ref_it) {
            std::string refName = (*ref_it).first;
            for (int useID = 0; useID < refToSameArrayCnt[refName]; useID++) {
                for (int reuseID = 0; reuseID < refToSameArrayCnt[refName]; reuseID++) {
                    
                    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops;
                    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops = findRefLoops(LoopRefTree, refName, useID, loops);
                    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops = findRefLoops(LoopRefTree, refName, reuseID, loops);
                    
                    std::string useIndvs = "";
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = useLoops.begin(), eit = useLoops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            useIndvs += "int ";
                            useIndvs += indvName[(*(*it)->LIS->IDV)[i]];
                            useIndvs += ", ";
                        }
                    }
                    useIndvs.pop_back();
                    useIndvs.pop_back();
                    
                    std::string reuseIndvs = "";
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = reuseLoops.begin(), eit = reuseLoops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            reuseIndvs += "int ";
                            reuseIndvs += indvName[(*(*it)->LIS->IDV)[i]];
                            reuseIndvs += "reuse";
                            reuseIndvs += ", ";
                        }
                    }
                    reuseIndvs.pop_back();
                    reuseIndvs.pop_back();
                    
                    errs() << "bool checkInterven";
                    errs() << refName + std::to_string(useID) + "_" + std::to_string(reuseID);
                    errs() << "(";
                    errs() << useIndvs + ", ";
                    errs() << reuseIndvs + ") { \n";
                    
                    checkIntervenBodyGen(LoopRefTree, refName, useID, reuseID);
                    
                    errs() << "    return false;\n";
                    errs() << "}\n";
                }
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::checkLocGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int useID, int reuseID) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops;
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops = findRefLoops(LoopRefTree, refName, useID, loops);
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops = findRefLoops(LoopRefTree, refName, reuseID, loops);
        
        std::string space = "    ";
        for (unsigned long i = 0; i < useLoops.size() + reuseLoops.size(); i++) {
            space += "    ";
        }
        
        /* Generate if statment to make sure they access the same location */
        errs() << space + "if ( calAddr" + refName + std::to_string(useID) + "(";
        std::string useIndvs = "";
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = useLoops.begin(), eit = useLoops.end(); it != eit; ++it) {
            for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                useIndvs += indvName[(*(*it)->LIS->IDV)[i]];
                useIndvs += ", ";
            }
        }
        useIndvs.pop_back();
        useIndvs.pop_back();
        
        errs() << useIndvs + ")";
        errs() << " == ";
        errs() << "calAddr" + refName + std::to_string(reuseID) + "(";
        
        std::string reuseIndvs = "";
        for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = reuseLoops.begin(), eit = reuseLoops.end(); it != eit; ++it) {
            for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                reuseIndvs += indvName[(*(*it)->LIS->IDV)[i]];
                reuseIndvs += "reuse";
                reuseIndvs += ", ";
            }
        }
        reuseIndvs.pop_back();
        reuseIndvs.pop_back();
        errs() << reuseIndvs + ") ) {\n";
        
        /* Generate chechInterven() function call */
        space += "    ";
        errs() << space +"if ( ";
        errs() << "checkInterven";
        errs() << refName + std::to_string(useID) + "_" + std::to_string(reuseID);
        errs() << "(" + useIndvs + ", " + reuseIndvs + ")";
        errs() << " == false) {\n";
        
        errs() << space + "  rtHistoCal(  rtCal" + refName + std::to_string(useID) + "_" + std::to_string(reuseID) + "(";
        errs() << useIndvs + ", " + reuseIndvs + ") );\n";
        
        errs() << space + "}\n";
        errs() << space + "findReuseFlag = true;\n";
        
        space.pop_back();
        space.pop_back();
        space.pop_back();
        space.pop_back();
        errs() << space + "}\n";
        
        /* finilize pair function */
        for (unsigned long i = 0; i < reuseLoops.size(); i++) {
            
            errs() << space + "if (findReuseFlag == true) {\n";
            errs() << space + "    break;\n";
            errs() << space + "}\n";
            
            space.pop_back();
            space.pop_back();
            space.pop_back();
            space.pop_back();
            errs() << space + "}\n";
        }
#if (SAMPLING == 0 || SAMPLING == 1)
        for (unsigned long i = 0; i < useLoops.size(); i++) {
            space.pop_back();
            space.pop_back();
            space.pop_back();
            space.pop_back();
            errs() << space + "}\n";
        }
#elif SAMPLING == 2
        space.pop_back();
        space.pop_back();
        space.pop_back();
        space.pop_back();
        errs() << space + "}\n";
#endif
        return;
    }
    
    void StaticSamplingCodeGen::pairRefRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int useID, int reuseID) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops;
        useLoopGen(LoopRefTree, refName, useID, loops);
        
        loops.clear();
        if (reuseLoopGen(LoopRefTree, refName, useID, reuseID, loops)) {
            checkLocGen(LoopRefTree, refName, useID, reuseID);
        }
        
        return;
    }
    
    // arrayName not used
    void StaticSamplingCodeGen::pairRefRTTopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        for (std::map<std::string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < it->second; i++) {
                for (int j = 0; j < it->second; j++) {
                    errs() << "void pair" + it->first + std::to_string(i) + "_" + std::to_string(j) + "() {\n";
                    pairRefRTBodyGen(LoopRefTree, it->first, i, j);
                    errs() << "}\n";
                }
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::rtCalFuncBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops, int useID, int reuseID) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator useIt = useLoops.begin();
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator reuseIt = reuseLoops.begin();
        
        // lowest common ancestor of use and reuse node
        int lcaIdx = -1;
        for (unsigned long i = 0; i < useLoops.size() && i < reuseLoops.size(); i++) {
            if (useLoops[i] == reuseLoops[i]) {
                if (i + 1 < useLoops.size() && i + i < reuseLoops.size()) {
                    if (useLoops[i+1] != reuseLoops[i+1]) {
                        lcaIdx = i;
                        break;
                    }
                } else {
                    lcaIdx = i;
                }
            }
        }
        
        std::string calExpr = "";
        
        while (useIt != useLoops.end() && reuseIt != reuseLoops.end()) {
            if ((*useIt) == (*reuseIt)) {
                calExpr += "(";
                calExpr += indvName[(*(*reuseIt)->LIS->IDV)[0]];
                calExpr += "reuse";
                calExpr += " - ";
                calExpr += indvName[(*(*useIt)->LIS->IDV)[0]];
                calExpr += ")";
                calExpr += " * ";
                calExpr += std::to_string(refCntOfLoop[*useIt]);
            } else {
                calExpr += "(";
                calExpr += getBound((*(*useIt)->LIS->LB)[0].second);
                calExpr += " - ";
                calExpr += indvName[(*(*useIt)->LIS->IDV)[0]];
                calExpr += ")";
                calExpr += " * ";
                calExpr += std::to_string(refCntOfLoop[*useIt]);
                
                calExpr += " + ";
                
                calExpr += "(";
                calExpr += indvName[(*(*reuseIt)->LIS->IDV)[0]];
                calExpr += "reuse";
                calExpr += " - ";
                calExpr += getBound((*(*reuseIt)->LIS->LB)[0].first);
                calExpr += ")";
                calExpr += " * ";
                calExpr += std::to_string(refCntOfLoop[*reuseIt]);
            }
            
            calExpr += " + ";
            
            ++useIt;
            ++reuseIt;
        }
        
        /* need to add the ref cnt between two different loops */
        //errs() << lcaIdx << " " << useLoops.size() << " " << reuseLoops.size() <<  "\n";
        if (lcaIdx != -1 && !(lcaIdx + 1 == (int) useLoops.size() && lcaIdx + 1 == (int) reuseLoops.size()) ) {
            
            //errs() << "here\n";
            
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *useRefNode = NULL;
            if (lcaIdx + 1 == (int) useLoops.size()) {
                useRefNode = findRef(LoopRefTree, refName, useID);
            }
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *reuseRefNode = NULL;
            if (lcaIdx + 1 == (int) reuseLoops.size()) {
                reuseRefNode = findRef(LoopRefTree, refName, reuseID);
            }
            bool useFlag = false;
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = useLoops[lcaIdx]->next->begin(), eit = useLoops[lcaIdx]->next->end(); it != eit; ++it) {
                if (reuseRefNode != NULL) {
                    if (reuseRefNode == (*it)) {
                        break;
                    }
                } else {
                    if (std::find(reuseLoops.begin(), reuseLoops.end(), *it) != reuseLoops.end()) {
                        break;
                    }
                }
                if (useFlag == true) {
                    if ((*it)->AA != NULL) {
                        calExpr += "1 + ";
                    }
                    if ((*it)->L != NULL) {
                        calExpr += std::to_string(refCntOfLoop[*it]);
                        calExpr += " + ";
                    }
                }
                if (useRefNode != NULL) {
                    if (useRefNode == (*it)) {
                        calExpr += " 1 + ";
                        useFlag = true;
                    }
                } else {
                    if (std::find(useLoops.begin(), useLoops.end(), *it) != useLoops.end()) {
                        useFlag = true;
                    }
                }
            }
        }
        
        if (lcaIdx == -1) {
            if (useLoops.size() == 0 && reuseLoops.size() != 0) {
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *useRefNode = findRef(LoopRefTree, refName, useID);
                bool useFlag = false;
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                    if (std::find(reuseLoops.begin(), reuseLoops.end(), *it) != reuseLoops.end()) {
                        
                        // calculate inside reuse loop
                        for (unsigned long i = 0; i + 1 < reuseLoops.size(); i++) {
                            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator reuseIter = reuseLoops[i]->next->begin(), reuseEIter = reuseLoops[i]->next->end(); reuseIter != reuseEIter; ++reuseIter) {
                                if ((*reuseIter) != reuseLoops[i+1]) {
                                    if ((*reuseIter)->AA != NULL) {
                                        calExpr += "1 + ";
                                    }
                                    if ((*reuseIter)->L != NULL) {
                                        calExpr += std::to_string(refCntOfLoop[*reuseIter]);
                                        calExpr += " * (";
                                        calExpr += getBound((*(*reuseIter)->LIS->LB)[0].second);
                                        calExpr += " - ";
                                        calExpr += getBound((*(*reuseIter)->LIS->LB)[0].first);
                                        calExpr += ")";
                                        calExpr += " + ";
                                    }
                                }
                            }
                        }
                        
                        break;
                    }
                    if (useFlag == true) {
                        if ((*it)->AA != NULL) {
                            calExpr += "1 + ";
                        }
                        if ((*it)->L != NULL) {
                            calExpr += std::to_string(refCntOfLoop[*it]);
                            calExpr += " + ";
                        }
                    }
                    if (useRefNode == (*it)) {
                        calExpr += " 1 + ";
                        useFlag = true;
                    }
                }
            }
            if (useLoops.size() != 0 && reuseLoops.size() == 0) {
                bool useFlag = false;
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *reuseRefNode = findRef(LoopRefTree, refName, reuseID);
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                    if (reuseRefNode == (*it)) {
                        break;
                    }
                    if (useFlag == true) {
                        if ((*it)->AA != NULL) {
                            calExpr += "1 + ";
                        }
                        if ((*it)->L != NULL) {
                            calExpr += std::to_string(refCntOfLoop[*it]);
                            calExpr += " + ";
                        }
                    }
                    if (std::find(useLoops.begin(), useLoops.end(), *it) != useLoops.end()) {
                        
                        // calculate inside use loop
                        
                        for (unsigned long i = 0; i + 1 < useLoops.size(); i++) {
                            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator useIter = useLoops[i]->next->begin(), useEIter = useLoops[i]->next->end(); useIter != useEIter; ++useIter) {
                                if ((*useIter) != useLoops[i+1]) {
                                    if ((*useIter)->AA != NULL) {
                                        calExpr += "1 + ";
                                    }
                                    if ((*useIter)->L != NULL) {
                                        calExpr += std::to_string(refCntOfLoop[*useIter]);
                                        calExpr += " * (";
                                        calExpr += getBound((*(*useIter)->LIS->LB)[0].second);
                                        calExpr += " - ";
                                        calExpr += getBound((*(*useIter)->LIS->LB)[0].first);
                                        calExpr += ")";
                                        calExpr += " + ";
                                    }
                                }
                            }
                        }
                        
                        useFlag = true;
                    }
                }
            }
            if (useLoops.size() != 0 && reuseLoops.size() != 0) {
                bool useFlag = false;
                
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                    if (std::find(reuseLoops.begin(), reuseLoops.end(), *it) != reuseLoops.end()) {
                        
                        // calculate inside reuse loop
                        for (unsigned long i = 0; i + 1 < reuseLoops.size(); i++) {
                            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator reuseIter = reuseLoops[i]->next->begin(), reuseEIter = reuseLoops[i]->next->end(); reuseIter != reuseEIter; ++reuseIter) {
                                if ((*reuseIter) != reuseLoops[i+1]) {
                                    if ((*reuseIter)->AA != NULL) {
                                        calExpr += "1 + ";
                                    }
                                    if ((*reuseIter)->L != NULL) {
                                        calExpr += std::to_string(refCntOfLoop[*reuseIter]);
                                        calExpr += " * (";
                                        calExpr += getBound((*(*reuseIter)->LIS->LB)[0].second);
                                        calExpr += " - ";
                                        calExpr += getBound((*(*reuseIter)->LIS->LB)[0].first);
                                        calExpr += ")";
                                        calExpr += " + ";
                                    }
                                }
                            }
                        }

                        
                        break;
                    }
                    if (useFlag == true) {
                        if ((*it)->AA != NULL) {
                            calExpr += "1 + ";
                        }
                        if ((*it)->L != NULL) {
                            calExpr += std::to_string(refCntOfLoop[*it]);
                            calExpr += " + ";
                        }
                    }
                    
                    if (std::find(useLoops.begin(), useLoops.end(), *it) != useLoops.end()) {
                        
                        // calculate inside use loop
                        for (unsigned long i = 0; i + 1 < useLoops.size(); i++) {
                            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator useIter = useLoops[i]->next->begin(), useEIter = useLoops[i]->next->end(); useIter != useEIter; ++useIter) {
                                if ((*useIter) != useLoops[i+1]) {
                                    if ((*useIter)->AA != NULL) {
                                        calExpr += "1 + ";
                                    }
                                    if ((*useIter)->L != NULL) {
                                        calExpr += std::to_string(refCntOfLoop[*useIter]);
                                        calExpr += " * (";
                                        calExpr += getBound((*(*useIter)->LIS->LB)[0].second);
                                        calExpr += " - ";
                                        calExpr += getBound((*(*useIter)->LIS->LB)[0].first);
                                        calExpr += ")";
                                        calExpr += " + ";
                                    }
                                }
                            }
                        }
                        
                        useFlag = true;
                    }
                }
            }
            if (useLoops.size() == 0 && reuseLoops.size() == 0) {
                bool useFlag = false;
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *useRefNode = findRef(LoopRefTree, refName, useID);
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *reuseRefNode = findRef(LoopRefTree, refName, reuseID);
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                    if (reuseRefNode == (*it)) {
                        break;
                    }
                    if (useFlag == true) {
                        if ((*it)->AA != NULL) {
                            calExpr += "1 + ";
                        }
                        if ((*it)->L != NULL) {
                            calExpr += std::to_string(refCntOfLoop[*it]);
                            calExpr += " + ";
                        }
                    }
                    if (useRefNode == (*it)) {
                        calExpr += " 1 + ";
                        useFlag = true;
                    }
                }
            }
        }
        
        /* to deal with the loop lengths are different */
        while (useIt != useLoops.end()) {
            calExpr += "(";
            calExpr += getBound((*(*useIt)->LIS->LB)[0].second);
            calExpr += " - ";
            calExpr += indvName[(*(*useIt)->LIS->IDV)[0]];
            calExpr += ")";
            calExpr += " * ";
            calExpr += std::to_string(refCntOfLoop[*useIt]);
            
            calExpr += " + ";
            
            useIt++;
        }
        
        while (reuseIt != reuseLoops.end()) {
            calExpr += "(";
            calExpr += indvName[(*(*reuseIt)->LIS->IDV)[0]];
            calExpr += "reuse";
            calExpr += " - ";
            calExpr += getBound((*(*reuseIt)->LIS->LB)[0].first);
            calExpr += ")";
            calExpr += " * ";
            calExpr += std::to_string(refCntOfLoop[*reuseIt]);
            
            calExpr += " + ";
            
            reuseIt++;
        }
        
        calExpr.pop_back();
        calExpr.pop_back();
        
        /* take ref access order into account */
        if ( !reuseLoops.empty() ) {
            if (reuseLoops.back() == useLoops.back() || std::find(useLoops.begin(), useLoops.end(), reuseLoops.back()) == useLoops.end() ) {
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = (*reuseLoops.rbegin())->next->begin(), eit = (*reuseLoops.rbegin())->next->end(); it != eit; ++it) {
                    if ((*it)->AA != NULL) {
                        if (refNumber[(*it)->AA] == reuseID && arrayName[(*it)->AA] == refName) {
                            calExpr += "+ ";
                            calExpr += std::to_string(refOrder[*it]);
                        }
                    }
                }
            }
        }
        
        if ( !useLoops.empty() ) {
            if (reuseLoops.back() == useLoops.back() || std::find(reuseLoops.begin(), reuseLoops.end(), useLoops.back()) == reuseLoops.end()) {
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = (*useLoops.rbegin())->next->begin(), eit = (*useLoops.rbegin())->next->end(); it != eit; ++it) {
                    if ((*it)->AA != NULL) {
                        if (refNumber[(*it)->AA] == useID && arrayName[(*it)->AA] == refName) {
                            calExpr += " - ";
                            calExpr += std::to_string(refOrder[*it]);
                        }
                    }
                }
            }
        }
        
        errs() << "    return " + calExpr + ";\n";
        
        return;
    }
    
    void StaticSamplingCodeGen::rtCalFuncTopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        for (std::map<std::string, int>::iterator ref_it = refToSameArrayCnt.begin(), ref_eit = refToSameArrayCnt.end(); ref_it != ref_eit; ++ref_it) {
            std::string refName = (*ref_it).first;
            for (int useID = 0; useID < refToSameArrayCnt[refName]; useID++) {
                for (int reuseID = 0; reuseID < refToSameArrayCnt[refName]; reuseID++) {
                    
                    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops;
                    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops = findRefLoops(LoopRefTree, refName, useID, loops);
                    std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops = findRefLoops(LoopRefTree, refName, reuseID, loops);
                    
                    std::string useIndvs = "";
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = useLoops.begin(), eit = useLoops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            useIndvs += "int ";
                            useIndvs += indvName[(*(*it)->LIS->IDV)[i]];
                            useIndvs += ", ";
                        }
                    }
                    useIndvs.pop_back();
                    useIndvs.pop_back();
                    
                    std::string reuseIndvs = "";
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = reuseLoops.begin(), eit = reuseLoops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            reuseIndvs += "int ";
                            reuseIndvs += indvName[(*(*it)->LIS->IDV)[i]];
                            reuseIndvs += "reuse";
                            reuseIndvs += ", ";
                        }
                    }
                    reuseIndvs.pop_back();
                    reuseIndvs.pop_back();
                    
                    errs() << "int rtCal" + refName + std::to_string(useID) + "_" + std::to_string(reuseID) + "(";
                    errs() << useIndvs + ", " + reuseIndvs;
                    errs() << ") {\n";
                    
                    rtCalFuncBodyGen(LoopRefTree, refName, useLoops, reuseLoops, useID, reuseID);
                    
                    errs() << "}\n";
                    
                }
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<std::string> indvs) {
        
        if (LoopRefTree->L != NULL) {
            for (std::vector<Value*>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvs.push_back(indvName[(*it)]);
            }
        }
        
        if (LoopRefTree->AA != NULL) {
            errs() << "int calAddr" + arrayName[LoopRefTree->AA] + std::to_string(refNumber[LoopRefTree->AA]) + "( ";
            std::string arguments = "";
            for (std::vector<std::string>::iterator it = indvs.begin(), eit = indvs.end(); it != eit; ++it) {
                arguments += "int " + (*it) + ", ";
            }
            arguments.pop_back();
            arguments.pop_back();
            errs() << arguments;
            errs() << ") {\n";
            errs() << "    int result = " + arrayExpression[LoopRefTree->AA] + ";\n";
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
    
    void StaticSamplingCodeGen::addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        std::vector<std::string> indvs;
        addrCalFuncGen(LoopRefTree, indvs);
        
        return;
    }
    
    void StaticSamplingCodeGen::headerGen() {
        
        errs() << "#include <map>\n";
        errs() << "#include <set>\n";
        errs() << "#include <cstdlib>\n";
        errs() << "#include <iostream>\n";
        errs() << "using namespace std;\n";
        
        return;
    }
    
    void StaticSamplingCodeGen::rtHistoGen() {
        
        errs() << "std::map<uint64_t, double> RT;\n";
        errs() << "std::map<uint64_t, double> MR;\n";
        errs() << "void rtHistoCal( int rt) {\n";
        errs() << "    if (RT.find(rt) == RT.end()) { \n";
        errs() << "        RT[rt] = 1;\n";
        errs() << "    } else {\n";
        errs() << "        RT[rt] += 1;\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        
        return;
    }
    
    void StaticSamplingCodeGen::rtDumpGen() {
        
        errs() << "void rtDump() {\n";
        errs() << "    cout << \"Start to dump reuse time histogram\\n\";\n";
        errs() << "    for (map<uint64_t, double>::iterator it = RT.begin(), eit = RT.end(); it != eit; ++it) {\n";
        errs() << "        cout << it->first << \" \" << it->second << \"\\n\";\n";
        errs() << "    }\n";
        errs() << "    return;\n";
        errs() << "}\n";
        
        return;
    }
    
    void StaticSamplingCodeGen::rtToMRGen() {
        
        errs() << "void RTtoMR_AET() {\n";
        
        std::string space = "    ";
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
    
    void StaticSamplingCodeGen::mrDumpGen() {
        
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
    
    void StaticSamplingCodeGen::mainGen() {
        
        errs() << "int main() {\n";
        
        std::string space = "    ";
        
        for (std::map<std::string, int>::iterator it = refToSameArrayCnt.begin(), eit = refToSameArrayCnt.end(); it != eit; ++it) {
            for (int i = 0; i < (*it).second; i++) {
                for (int j = 0; j < (*it).second; j++) {
                    //errs() << space + "cout << \" check pair " + (*it).first + std::to_string(i) + " " + (*it).first + std::to_string(j) + "\\n \";\n";
                    errs() << space + "pair" + (*it).first + std::to_string(i) + "_" + std::to_string(j) + "();\n";
                    //errs() << "    rtDump();\n";
                }
            }
        }
        
        //errs() << "    rtDump();\n";
        
        errs() << "    RTtoMR_AET();";
        
        errs() << "    dumpMR();";
        
        errs() << "    return 0;\n";
        errs() << "}\n";
        
        return;
    }
    
    void StaticSamplingCodeGen::numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
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
                numberRefToSameArray( (*it) );
            }
        }
        
        return;
    }
    
    int StaticSamplingCodeGen::initRefTotalOrder(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, int order) {
        
        if (LoopRefTree == NULL) {
            return -1;
        }
        
        if (LoopRefTree->AA != NULL) {
            refTotalOrder[LoopRefTree] = order;
            
            //errs() << refTotalOrder[LoopRefTree] << " " << refOrder[LoopRefTree] << "\n";
            
            order++;
            return order;
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                order = initRefTotalOrder(*it, order);
            }
        }
        
        return order;
    }
    
    int StaticSamplingCodeGen::initRefCntOfLoop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree == NULL) {
            return 0;
        }
        
        if (LoopRefTree->L != NULL && LoopRefTree->next != NULL) {
            
            int cnt = 0;
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                
                if ((*it)->AA != NULL) {
                    cnt += 1;
                }
                
                if ((*it)->L != NULL) {
                    
                    /* need to handle symbolic bounds */
                    
                    int lB = stoi(getBound((*(*it)->LIS->LB)[0].first));
                    int uB = stoi(getBound((*(*it)->LIS->LB)[0].second));
                    cnt += initRefCntOfLoop(*it) * (uB - lB);
                }
            }
            
            refCntOfLoop[LoopRefTree] = cnt;
            return cnt;
            
        }
        
        /* root node */
        if (LoopRefTree->L == NULL && LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                if ((*it)->L != NULL) {
                    initRefCntOfLoop(*it);
                }
            }
        }
        
        return 0;
    }
    
    void StaticSamplingCodeGen::initLoopOrder(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree == NULL) {
            return;
        }
        
        if (LoopRefTree->next != NULL) {
            int order = 0;
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                if ((*it)->L != NULL) {
                    loopOrder[*it] = order;
                    order++;
                    initLoopOrder(*it);
                }
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::initRefOrder(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
        if (LoopRefTree == NULL) {
            return;
        }
        
        if (LoopRefTree->next != NULL) {
            int order = 0;
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                if ((*it)->AA != NULL) {
                    refOrder[*it] = order;
                    order++;
                }
                if ((*it)->L != NULL) {
                    initRefOrder(*it);;
                }
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree) {
        
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
    
    void StaticSamplingCodeGen::initArrayName() {
        
        for (std::map<Instruction*, std::string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            it->second.replace(std::find(it->second.begin(), it->second.end(), '.'), std::find(it->second.begin(), it->second.end(), '.') +1, 1, '_');
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::init(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {
        
        initLoopOrder(LoopRefTree);
        
        initRefOrder(LoopRefTree);
        
        initRefCntOfLoop(LoopRefTree);
        
        initIndvName(LoopRefTree);
        
        initArrayName();
        
        initRefTotalOrder(LoopRefTree, 0);
        
        return;
    }
    
    bool StaticSamplingCodeGen::runOnFunction(Function &F) {
        
        errs() << " // Start to generating Static Sampling Code\n";
        
        /* reading info from previous passes */
        arrayName = getAnalysis<idxAnalysis::IndexAnalysis>().arrayName;
        arrayExpression = getAnalysis<idxAnalysis::IndexAnalysis>().arrayExpression;
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree;
        
        
        if (LoopRefTree == NULL) {
            errs() << "No array access\n";
            return false;
        }
        
        init(LoopRefTree);
        
        /* number the references to the same array */
        numberRefToSameArray(LoopRefTree);
        
        /* generate headers */
        headerGen();
        
        /* generate rtHistoCal function */
        rtHistoGen();
        
        /* generate rtToMR function */
        rtToMRGen();
        
        /* generate address calculation functions */
        addrCalFuncGenTop(LoopRefTree);
        
        /* generate reuse time calculation function */
        rtCalFuncTopGen(LoopRefTree);
        
        /* generate checkInterven top function */
        checkIntervenTopGen(LoopRefTree);
        
        /* generate pair-wise rt top function */
        pairRefRTTopGen(LoopRefTree);
        
        /* generate rtDump function */
        rtDumpGen();
        
        /* generate mrDump function */
        mrDumpGen();
        
        /* generate main function */
        mainGen();
        
        return false;
    }
    
    void StaticSamplingCodeGen::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
        
        return;
    }
}
