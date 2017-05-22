//
//  gVarAnalysis.hpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#ifndef gVarAnalysis_hpp
#define gVarAnalysis_hpp

#include <stdio.h>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Analysis/CFG.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"

using namespace llvm;

namespace gVarAnalysis {
    struct GlobalVariableAnalysis : public FunctionPass {
        static char ID;
        GlobalVariableAnalysis();
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
        
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}



#endif /* gVarAnalysis_hpp */
