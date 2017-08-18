//
//  ssCodeGen.hpp
//  LLVM
//
//  Created by dongchen on 8/16/17.
//
//

#ifndef ssCodeGen_hpp
#define ssCodeGen_hpp

#include <map>

#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
//#include "brchAnalysis.hpp"

using namespace llvm;

namespace ssCodeGen {
    
    struct StaticSamplingCodeGen : public FunctionPass {
        
        static char ID;
        StaticSamplingCodeGen();
        
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
        
        std::map<std::string, int> refToSameArrayCnt;
        std::map<Instruction*, int> refNumber;
        void numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        void addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<std::string> indvs);
        void addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        void pairRefRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int useID, int reuseID);
        void pairRefRTTopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        void useLoopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> indvs);
        std::string getBound(Value *bound);
        
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
        
    };
}

#endif /* ssCodeGen_hpp */
