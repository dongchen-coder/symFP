//
//  loopAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#include "loopAnalysis.hpp"

using namespace llvm;
using namespace std;

namespace loopAnalysis {
    
    char LoopIndvBoundAnalysis::ID = 0;
    static RegisterPass<LoopIndvBoundAnalysis> X("loopAnalysis", "loop induction variable/bound analysis Pass", false, false);
    
    LoopIndvBoundAnalysis::LoopIndvBoundAnalysis() : FunctionPass(ID) {}
    
    
    
    LoopIndvBoundAnalysis::LoopBound LoopIndvBoundAnalysis::findLoopBound(Loop *L, Value *var) {
        Value* ub;
        Value* lb;
        for (BasicBlock *b : L->getBlocks()) {
            if (std::regex_match(b->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)"))) {
                for(BasicBlock::iterator it = b->begin(), eit = b->end(); it != eit; ++it) {
                    if (isa<ICmpInst>(*it)) {
//                        errs() << "Upper Bound for Loop is: ";
                        Value *v = it->getOperand(1);
                        if (isa<ConstantInt>(v)) {
//                            errs() << dyn_cast<ConstantInt>(v)->getValue();
                            ub = v;
                        }
                    }
                }
                // iterate all its current predecessors
                for(auto pit = pred_begin(b), pet = pred_end(b); pit != pet; ++pit) {
                    BasicBlock *pred = *pit;
                    if (!std::regex_match(pred->getName().str(), std::regex("^for.inc\\d*$|^for.inc$"))) {
                        // iterate all store instructions, find the first operand when the second operand is a idv
                        for (BasicBlock::iterator it = pred->begin(), eit = pred->end(); it != eit; ++it) {
                            if (isa<StoreInst>(*it)) {
                                Value *v = it->getOperand(0);
                                if (it->getOperand(1)->getName().equals(var->getName())) {
//                                    errs() << "Lower Bound for Loop is: " << dyn_cast<ConstantInt>(v)->getValue() <<"\n";
                                    lb = v;
                                }
                            }
                        }
                    }
                }
            }
            
        }
        return make_pair(lb, ub);
    }

    string LoopIndvBoundAnalysis::getBound(Value *bound) {
        
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

    
    
    
    /*
     * 获得Loop L所有子Loop中的cond BasicBlock
     */
    
    vector<BasicBlock *> LoopIndvBoundAnalysis::getSubLoopCondBlock(Loop *L) {
        vector<BasicBlock *> temp;
        for (Loop *SL : L->getSubLoops()) {
            for (BasicBlock *BB: SL->getBlocks()) {
                if (std::regex_match (BB->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)"))) {
                    //                    errs() << BB->getName() + " \n";
                    temp.push_back(BB);
                }
            }
        }
        return temp;
    }
    
    
    LoopIndvBoundAnalysis::LoopInfoStruct* LoopIndvBoundAnalysis::ExtractLoopInfo(Loop *L) {
        
        vector<BasicBlock *> subLoopCondBlocks = getSubLoopCondBlock(L);
        
        vector<Value *>* IDV = new vector<Value *>;
        vector<LoopIndvBoundAnalysis::LoopBound>* LB = new vector<LoopIndvBoundAnalysis::LoopBound>;
        
        LoopIndvBoundAnalysis::LoopInfoStruct* LIS = new LoopIndvBoundAnalysis::LoopInfoStruct;
        
        // find all BasicBlocks in the loop L
        for (BasicBlock *b : L->getBlocks()) {
            if (std::regex_match (b->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)")) && find(subLoopCondBlocks.begin(), subLoopCondBlocks.end(), b) == subLoopCondBlocks.end()) {
                for (Instruction &II : *b) {
                    if (isa<LoadInst>(II)) {
                        IDV->push_back(II.getOperand(0));
                        LB->push_back(findLoopBound(L, II.getOperand(0)));
                    }
                }
            }
        }
        
        LIS->IDV = IDV;
        LIS->LB = LB;
        
        return LIS;
    }
    
    void LoopIndvBoundAnalysis::LoopTreeConstruction(Loop* L, LoopTreeNodes * root, int level) {
        
        if (L != NULL) {
            
            //root = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
            root->LoopLevel = level;
            root->LIS = ExtractLoopInfo(L);
            root->next = new vector<LoopTreeNodes *>;
            
            for (Loop *SL : L->getSubLoops()) {
                LoopTreeNodes * subTmp = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
                root->next->push_back(subTmp);
                LoopTreeConstruction(SL, subTmp, level+1);
            }
        }
        
        return;
    }
    
    
    LoopIndvBoundAnalysis::LoopTreeNodes* LoopIndvBoundAnalysis::LoopTreeConstructionTop(LoopTreeNodes* root) {
        
        LoopInfo &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
        
        if (!LI.empty()) {
            root = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
            root->LoopLevel = 0;
            root->LIS = NULL;
            root->next = new vector<LoopTreeNodes *>;
            
            for(LoopInfo::iterator it = LI.begin(), eit = LI.end(); it != eit; ++it){
                
                LoopTreeNodes * Tmp = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
                Tmp->LIS = ExtractLoopInfo(*it);
                Tmp->LoopLevel = 1;
                Tmp->next = new vector<LoopTreeNodes *>;
                
                root->next->push_back(Tmp);
                
                for (Loop *SL : (*it)->getSubLoops()) {
                    LoopTreeNodes * subTmp = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
                    Tmp->next->push_back(subTmp);
                    LoopTreeConstruction(SL, subTmp, 2);
                }
            }
        }
        
        return root;
    }
    
    void LoopIndvBoundAnalysis::DumpLoopTree(LoopTreeNodes* LTroot, std::string prefix) {
        
        if (LTroot != NULL) {
            
            if (LTroot->LIS != NULL) {
                for (vector<Value*>::iterator it = LTroot->LIS->IDV->begin(), eit = LTroot->LIS->IDV->end(); it != eit; ++it ) {
                    if ((*it) != NULL) {
                        errs() << prefix << (*it)->getName() << "\n";
                    } else {
                        errs() << prefix << "NULL idv\n";
                    }
                }
                
                for (vector<LoopBound>::iterator it = LTroot->LIS->LB->begin(), eit = LTroot->LIS->LB->end(); it != eit; ++it) {
                    errs() << prefix << "Loop Bound: (" << getBound((*it).first) << ", " << getBound((*it).second) << ")\n";
                }
            }
            
            if (LTroot->next != NULL) {
                for (vector<LoopTreeNodes*>::iterator it = LTroot->next->begin(), eit = LTroot->next->end(); it != eit; ++it) {
                    DumpLoopTree(*it, prefix + "--");
                }
            }
        } else {
            errs() << "LTroot == NULL\n";
        }
        
        return;
    }
    
    /* Main function for loop induction variable analysis and loop bound analysis */
    bool LoopIndvBoundAnalysis::runOnFunction(Function &F) {
        
        errs() << "\nStart analysis loops\n";
        
        LoopTreeNodes* LTroot = NULL;
        LTroot = LoopTreeConstructionTop(LTroot);
        
        DumpLoopTree(LTroot, "");
        
        errs() << "\nFinish analysis loops\n";
        
        return false;
    }
    
    /* Dependence relation of the analysis paths */
    void LoopIndvBoundAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        return;
    }
    
}
