//
//  gVarAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#include "gVarAnalysis.hpp"

namespace gVarAnalysis {
    char GlobalVariableAnalysis::ID = 0;
    static RegisterPass<GlobalVariableAnalysis> X("gVarAnalysis", "global variable analysis Pass", false, false);

    GlobalVariableAnalysis::GlobalVariableAnalysis() : FunctionPass(ID) {}
    
    bool GlobalVariableAnalysis::runOnFunction(Function &F) {
        
        errs() << "\n /* Start to analysis global variable \n";
        
        errs() << "\n Finish to analysis global variable */ \n";
        
        return false;
    }
    
    void GlobalVariableAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        return;
    }
}
