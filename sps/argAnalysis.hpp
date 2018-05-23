//
//  argAnalysis.hpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#ifndef argAnalysis_hpp
#define argAnalysis_hpp

#include <stdio.h>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Analysis/CFG.h"
#include "idxAnalysis.hpp"

using namespace llvm;

namespace argAnalysis {
    struct ArgumentAnalysis : public FunctionPass {
        static char ID;
        ArgumentAnalysis();
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
        
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif /* argAnalysis_hpp */
