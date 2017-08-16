//
//  loopAnalysis.hpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#ifndef loopAnalysis_hpp
#define loopAnalysis_hpp

#include <stdio.h>
#include <set>
#include <regex>
#include <map>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/CFG.h"
#include "llvm/IR/Constants.h"


#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"

using namespace llvm;
using namespace std;

#define LOOP_DEBUG

namespace loopAnalysis {
    struct LoopIndvBoundAnalysis : public FunctionPass {
        static char ID;
        LoopIndvBoundAnalysis();
        
        typedef pair<Value*, Value*> LoopBound;
        
        struct LoopInfoStruct{
            vector<Value *>* IDV;
            vector<LoopIndvBoundAnalysis::LoopBound>* LB;
        };
        
        /* The nodes in Loop tree can either be a loop or an array access */
        struct LoopTreeNodes {
            int LoopLevel;
            Loop* L;
            LoopInfoStruct* LIS;
            Instruction* AA;
            vector<LoopTreeNodes *>* next;
        };
        
        LoopTreeNodes* LoopInfoTree;
        
        //vector<LoopInfoStruct> LoopInfoVector;
        
        //void subLoop(Loop *L);
        /* Find all Basic Induction Variable */
        //void findIDV(Loop *L);
        
        /* Find the Loop Bound */
        LoopBound findLoopBound(Loop *L, Value *var);
        
        /* Find All the BasicBlocks within the SubLoops */
        vector<BasicBlock *> getSubLoopCondBlock(Loop *L);
        
        /* Print Func */
        //void dumpLoopInfoStruct();
        
        /* Get the Loop Bound in a string format */
        string getBound(Value *);
        
        void DumpLoopTree(LoopTreeNodes* LTroot, std::string prefix);
        LoopInfoStruct* ExtractLoopInfo(Loop *L);
        void LoopTreeConstruction(Loop* LI , LoopTreeNodes * root, int level);
        LoopTreeNodes* LoopTreeConstructionLoop(LoopTreeNodes * root);
        LoopTreeNodes* ExtractRefInfo(Instruction* Inst);
        LoopTreeNodes* LoopTreeConstructionRef(LoopTreeNodes * root, vector<BasicBlock*> BBList);
        
        //void inductionVariableAnalysis(Function &F);
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
        
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}



#endif /* loopAnalysis_hpp */
