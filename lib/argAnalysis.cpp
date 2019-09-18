//
//  argAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#include "argAnalysis.hpp"

using namespace llvm;

namespace argAnalysis {
    char ArgumentAnalysis::ID = 0;
    static RegisterPass<ArgumentAnalysis> X("argAnalysis", "Argument analysis Pass", false, false);
    
    ArgumentAnalysis::ArgumentAnalysis() : FunctionPass(ID) {}
    
    bool ArgumentAnalysis::runOnFunction(Function &F) {
        
        errs() << "\n /* Start to analyze argument\n";
        
        for (Function::arg_iterator it = F.arg_begin(), eit = F.arg_end(); it != eit; ++it) {
            it->dump();
        }
        
        errs() << "\n Start to analysis argument */ \n";
        
        return false;
    }
    
    void ArgumentAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        return;
    }
}
