//
//  loopAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#include <regex>
#include "loopAnalysis.hpp"

using namespace llvm;
using namespace std;

namespace loopAnalysis {
    
    char LoopIndvBoundAnalysis::ID = 0;
    static RegisterPass<LoopIndvBoundAnalysis> X("loopAnalysis", "loop induction variable/bound analysis Pass", false, false);
    
    
    struct LoopInfoStruct{
        Loop *L;
        Value * IDV;
        LoopIndvBoundAnalysis::LoopBound LB;
        vector<Loop *> SL;
    };
    
    vector<LoopInfoStruct> LoopInfoVector;
    
    LoopIndvBoundAnalysis::LoopIndvBoundAnalysis() : FunctionPass(ID) {}
    void LoopIndvBoundAnalysis::subLoop(Loop *L) {
        for (Loop *SL : L->getSubLoops()) {
            errs() << "Sub Loop:\n";
//            findIDV(SL);
            subLoop(SL);
        }
        
        return;
    }
    
    void LoopIndvBoundAnalysis::findIDV(Loop *L) {
        vector<BasicBlock *> subLoopCondBlocks = getSubLoopCondBlock(L);
        
        // find all BasicBlocks in the loop L
        for (BasicBlock *b : L->getBlocks()) {
            if (std::regex_match (b->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)")) && find(subLoopCondBlocks.begin(), subLoopCondBlocks.end(), b) == subLoopCondBlocks.end()) {
                for (Instruction &II : *b) {
                    if (isa<LoadInst>(II)) {
                        LoopInfoStruct temp = {
                            L,
                            II.getOperand(0),
                            findLoopBound(L, II.getOperand(0)),
                            L->getSubLoops()
                        };
                        
                        LoopInfoVector.push_back(temp);
                    }
                }
            }
        }
        
        errs() << "\n";
    }
    
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
                                if (isa<ConstantInt>(v) && it->getOperand(1)->getName().equals(var->getName())) {
//                                    errs() << "Lower Bound for Loop is: " << dyn_cast<ConstantInt>(v)->getValue() <<"\n";
                                    lb = v;
                                }
                                else {
                                    // TO DO
                                }
                            }
                        }

                    }
                }
            }
            
        }
        return make_pair(lb, ub);
    }
    
    void LoopIndvBoundAnalysis::inductionVariableAnalysis(Function &F) {
        errs() << "Induction variable analysis\n";
        
        LoopInfo &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
        if (!LI.empty()) {
            for(LoopInfo::iterator it = LI.begin(), eit = LI.end(); it != eit; ++it){
                errs() << "Loop:\n";
//                findIDV(*it);
                subLoop(*it);
            }
        }
        dumpLoopInfoStruct();
        return;
    }
    
    bool LoopIndvBoundAnalysis::runOnFunction(Function &F) {
        
        errs() << "\nStart analysis loops\n";
        
        inductionVariableAnalysis(F);
        
        return false;
    }
    
    void LoopIndvBoundAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        return;
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
    
    void LoopIndvBoundAnalysis::dumpLoopInfoStruct() {
        for (LoopInfoStruct is: LoopInfoVector) {
            errs() << "For Loop " << is.L->getName() << ": \n";
            errs() << "Induction Variable: " << is.IDV->getName() << "\n";
            errs() << "Loop Bound: (" << dyn_cast<ConstantInt>(is.LB.first)->getValue() << ", " << dyn_cast<ConstantInt>(is.LB.second)->getValue() << ")\n";
            errs() << "Sub Loops: ";
            for (Loop* sl: is.SL) {
                errs() << sl->getName() << " ";
            }
            errs() << "\n\n";
        }
    }
}
