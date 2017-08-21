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
                std::string space = "    ";
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops.begin(), elit = loops.end(); lit != elit; ++lit) {
                    for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                        errs() << space + "for ( int " + (*(*lit)->LIS->IDV)[i]->getName() << " = ";
                        errs() << getBound((*(*lit)->LIS->LB)[i].first) << "; ";
                        errs() << (*(*lit)->LIS->IDV)[i]->getName() << " < " + getBound((*(*lit)->LIS->LB)[i].second) + "; ";
                        errs() << (*(*lit)->LIS->IDV)[i]->getName() << "++) {\n";
                    }
                    space += "    ";
                }
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
    
    void StaticSamplingCodeGen::reuseLoopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int useID, int reuseID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops = findRefLoops(LoopRefTree, refName, useID, loops);
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops = findRefLoops(LoopRefTree, refName, reuseID, loops);
        
        bool breakflag = false;
        if (useLoops.size() != 0 && reuseLoops.size() != 0 && loopOrder[useLoops[0]] > loopOrder[reuseLoops[0]]) {
            breakflag = true;
        }
        
        std::string space = "    ";
        for (unsigned long i = 0; i < useLoops.size(); i++) {
            space += "    ";
        }
        
        if (breakflag == true) {
            for (unsigned long i = 0; i < useLoops.size(); i++) {
                space += "    ";
                errs() << space + "break;\n";
            }
            return;
        }
        
        for (unsigned long i = 0; i < reuseLoops.size(); i++) {
            
            if (loopOrder[useLoops[i]] == loopOrder[reuseLoops[i]]) {
                for (unsigned long j = 0; j < reuseLoops[i]->LIS->IDV->size(); j++) {
                    errs() << space + "for ( int " + (*reuseLoops[i]->LIS->IDV)[j]->getName() + "reuse";
                    errs() << " = ";
                    errs() << (*useLoops[i]->LIS->IDV)[j]->getName() << "; ";
                    errs() << (*reuseLoops[i]->LIS->IDV)[j]->getName() + "reuse";
                    errs() << " < ";
                    errs() << getBound((*reuseLoops[i]->LIS->LB)[j].second);
                    errs() << "; ";
                    errs() << (*reuseLoops[i]->LIS->IDV)[j]->getName() + "reuse" + "++";
                    errs() << ") {\n";
                    space += "    ";
                }
            } else {
                for (unsigned long j = 0; j < reuseLoops[i]->LIS->IDV->size(); j++) {
                    errs() << space + "for ( int " + (*reuseLoops[i]->LIS->IDV)[j]->getName() + "reuse";
                    errs() << " = ";
                    errs() << getBound((*reuseLoops[i]->LIS->LB)[j].first) << "; ";
                    errs() << (*reuseLoops[i]->LIS->IDV)[j]->getName() + "reuse";
                    errs() << " < ";
                    errs() << getBound((*reuseLoops[i]->LIS->LB)[j].second);
                    errs() << "; ";
                    errs() << (*reuseLoops[i]->LIS->IDV)[j]->getName() + "reuse" + "++";
                    errs() << ") {\n";
                    space += "    ";
                }
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::checkIntervenBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int useID, int reuseID) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops;
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops = findRefLoops(LoopRefTree, refName, useID, loops);
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops = findRefLoops(LoopRefTree, refName, reuseID, loops);
        
        
        
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
                            useIndvs += (*(*it)->LIS->IDV)[i]->getName();
                            useIndvs += ", ";
                        }
                    }
                    useIndvs.pop_back();
                    useIndvs.pop_back();
                    
                    std::string reuseIndvs = "";
                    for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = reuseLoops.begin(), eit = reuseLoops.end(); it != eit; ++it) {
                        for (unsigned long i = 0; i < (*it)->LIS->IDV->size(); i++) {
                            reuseIndvs += "int ";
                            reuseIndvs += (*(*it)->LIS->IDV)[i]->getName();
                            reuseIndvs += "reuse";
                            reuseIndvs += ", ";
                        }
                    }
                    reuseIndvs.pop_back();
                    reuseIndvs.pop_back();
                    
                    errs() << "void checkInterven";
                    errs() << refName + std::to_string(useID) + "_" + std::to_string(reuseID);
                    errs() << "(";
                    errs() << useIndvs + ", ";
                    errs() << reuseIndvs + ") { \n";
                    
                    checkIntervenBodyGen(LoopRefTree, refName, useID, reuseID);
                    
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
                useIndvs += (*(*it)->LIS->IDV)[i]->getName();
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
                reuseIndvs += (*(*it)->LIS->IDV)[i]->getName();
                reuseIndvs += "reuse";
                reuseIndvs += ", ";
            }
        }
        reuseIndvs.pop_back();
        reuseIndvs.pop_back();
        errs() << reuseIndvs + ") {\n";
        
        /* Generate chechInterven() function call */
        space += "    ";
        errs() << space + "checkInterven";
        errs() << refName + std::to_string(useID) + "_" + std::to_string(reuseID);
        errs() << "(" + useIndvs + ", " + reuseIndvs + ");\n";
        errs() << space + "break;\n";
        
        /* finilize pair function */
        for (unsigned long i = 0; i < useLoops.size() + reuseLoops.size() + 1; i++) {
            space.pop_back();
            space.pop_back();
            space.pop_back();
            space.pop_back();
            errs() << space + "}\n";
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::pairRefRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int useID, int reuseID) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops;
        useLoopGen(LoopRefTree, refName, useID, loops);
        
        loops.clear();
        reuseLoopGen(LoopRefTree, refName, useID, reuseID, loops);
        
        checkLocGen(LoopRefTree, refName, useID, reuseID);
        
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
    
    void StaticSamplingCodeGen::addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<std::string> indvs) {
        
        if (LoopRefTree->L != NULL) {
            for (std::vector<Value*>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
                indvs.push_back((*it)->getName());
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
                    initRefOrder(*it);
                }
            }
        }
        
        return;
    }
    
    void StaticSamplingCodeGen::init(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {
        
        initLoopOrder(LoopRefTree);
        
        initRefOrder(LoopRefTree);
        
        return;
    }
    
    bool StaticSamplingCodeGen::runOnFunction(Function &F) {
        
        errs() << "Start to generating Static Sampling Code\n";
        
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
        
        /* generate address calculation functions */
        addrCalFuncGenTop(LoopRefTree);
        
        /* generate pair-wise rt top function */
        pairRefRTTopGen(LoopRefTree);
        
        /* generate checkInterven top function */
        checkIntervenTopGen(LoopRefTree);
        
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
