//
//  loopAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#include "loopAnalysis.hpp"

using namespace llvm;

namespace loopAnalysis {
    
    char LoopIndvBoundAnalysis::ID = 0;
    static RegisterPass<LoopIndvBoundAnalysis> X("loopAnalysis", "loop induction variable/bound analysis Pass", false, false);
    
    LoopIndvBoundAnalysis::LoopIndvBoundAnalysis() : FunctionPass(ID) {}
    
    void LoopIndvBoundAnalysis::subLoop(Loop *L) {
        
        for (Loop *SL : L->getSubLoops()) {
            errs() << "Loop:\n";
            PHINode* phi = SL->getCanonicalInductionVariable();
            if (phi != NULL) {
                phi->dump();
            } else {
                errs() << "  Can not find induction variable\n";
            }
            subLoop(SL);
        }
        
        return;
    }
    
    void LoopIndvBoundAnalysis::inductionVariableAnalysis(Function &F) {
        errs() << "Induction variable analysis\n";
        
        LoopInfo &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
        
        if (!LI.empty()) {
            for(LoopInfo::iterator it = LI.begin(), eit = LI.end(); it != eit; ++it){
                errs() << "Loop:\n";
                PHINode* phi = (*it)->getCanonicalInductionVariable();
                if (phi != NULL) {
                    phi->dump();
                } else {
                    errs() << " Can not find induction variable\n";
                }
                subLoop(*it);
            }
        }
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
}
