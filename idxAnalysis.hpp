//
//  idxAnalysis.hpp
//  LLVM
//
//  Created by dongchen on 5/18/17.
//
//

#ifndef idxAnalysis_hpp
#define idxAnalysis_hpp

#include <stdio.h>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/Instructions.h"
#include <map>
#include <vector>

//#define IDX_DEBUG

using namespace llvm;

namespace idxAnalysis {
    struct IndexAnalysis : public FunctionPass {
        static char ID;
        IndexAnalysis();
        
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
        
        std::vector<Instruction*> instStackInit(Instruction* inst);
        std::string computeExpression(Instruction* inst);
        
        std::string getArrayName(GetElementPtrInst* inst);
        
        void findAllArrayAccesses(Function &F);
        
        void dumpAllInfo();
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
    };
}


#endif /* idxAnalysis_hpp */
