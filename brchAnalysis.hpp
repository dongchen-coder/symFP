//
//  brchAnalysis.hpp
//  LLVM
//
//  Created by Fangzhou Liu on 2017/7/11.
//
//

#ifndef brchAnalysis_hpp
#define brchAnalysis_hpp

#include <stdio.h>
#include <regex>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Analysis/CFG.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"

using namespace llvm;

namespace brchAnalysis {
    
    struct BranchAnalysis : public FunctionPass {
        static char ID;
        BranchAnalysis();
        
        /* Extract the Conditions From the Block Within the Loop */
        void BranchConditionAnalysis(Function &F);
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
        
        /* Find the If Body in the Loop */
        vector<BasicBlock *> FindIfBody(Loop *L);
        
        /* Find the Branch Conditions */
        void FindBranch(Loop *L);
        
        /* Check whether the branch is valid */
        bool CheckBranchValid(Loop *L, Value *v);

        /* Check Whether the Condition Var belongs to the self/parent Loop Indv */
        bool CheckIndvVar(Loop *L, string var);
        
        /* Print the Conditions */
        void dumpBranchInfoStruct();
        
        /* Return the String Representation of a Value */
        string dumpValue(Value *v);
        
        /* Return the String Representation of a Array */
        string dumpArray(Value *arr);
        
        
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}


#endif /* brchAnalysis_hpp */
