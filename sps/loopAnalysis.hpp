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

/* Debugging flag, enable if debug is needed */
//#define LOOP_DEBUG

namespace loopAnalysis {
    struct LoopIndvBoundAnalysis : public FunctionPass {
        static char ID;
        LoopIndvBoundAnalysis();
        
        /* array index info: init from idxAnalysis pass */
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
        
        /* Data structures for loop info */
        typedef pair<Value*, Value*> LoopBound;              /* loop bounds are stored in pairs (lower bound, upper bound) */
        
        struct LoopInfoStruct{
            vector<Value *>* IDV;                            /* induction variable */
            vector<Value *>* INC;                            /* stride */
            vector<LoopIndvBoundAnalysis::LoopBound>* LB;    /* bounds */
            vector<llvm::CmpInst::Predicate>* PREDICATE;     /* upper bound condition >, >=, <, <= */
        };
        
        struct LoopRefTNode {                                /* Tree node structure for loop/reference tree */
            int LoopLevel;
            Loop* L;
            LoopInfoStruct* LIS;
            Instruction* AA;
            vector<LoopRefTNode *>* next;
        };
        
        LoopRefTNode* LoopRefTree;                           /* Root node for loop reference tree */
        
        /* Find the Loop Bound */
        LoopBound findLoopBound(Loop *L, Value *var);
        
        /* Find All the Condition/Increment BasicBlocks within the SubLoops */
        vector<BasicBlock *> getSubLoopCondBlock(Loop *L);
        vector<BasicBlock *> getSubLoopIncBlock(Loop *L);
        
        /* Get the Loop Bound in a string format */
        string getBound(Value * expr);
        
        /* Converting predicate to string */
        string predicateToString(llvm::CmpInst::Predicate PREDICATE);
        
        /* Create and init node for references */
        LoopRefTNode* ExtractRefInfo(Instruction* Inst);
        
        /* Decroate the loop tree with references */
        LoopRefTNode* LoopTreeConstructionRef(LoopRefTNode * root, vector<BasicBlock*> BBList);
        
        /* Init loopInfoStruct for each loop L */
        LoopInfoStruct* ExtractLoopInfo(Loop *L);
        
        /* Construct loop tree (references are not included) */
        void LoopTreeConstruction(Loop* LI , LoopRefTNode * root, int level);     /* called by LoopTreeConstructionLoop() */
        LoopRefTNode* LoopTreeConstructionLoop(LoopRefTNode * root);              /* Top */
        
        /* Dump loop/reference tree */
        void DumpLoopTree(LoopRefTNode* LTroot, std::string prefix);
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;

        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}



#endif /* loopAnalysis_hpp */
