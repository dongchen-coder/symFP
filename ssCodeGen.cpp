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
    
    void StaticSamplingCodeGen::pairRefRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string arrayName, int useID, int reuseID) {
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops;
        useLoopGen(LoopRefTree, arrayName, useID, loops);
        
        loops.clear();
        //reuseLoopGen(LoopRefTree, arrayName, reuseID, loops);
        
        
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
        
        /* number the references to the same array */
        numberRefToSameArray(LoopRefTree);
        
        /* generate address calculation functions */
        addrCalFuncGenTop(LoopRefTree);
        
        /* generate pair-wise rt top function */
        pairRefRTTopGen(LoopRefTree);
        
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
