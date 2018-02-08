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

//#define LOOP_DEBUG

namespace loopAnalysis {
    struct LoopIndvBoundAnalysis : public FunctionPass {
        static char ID;
        LoopIndvBoundAnalysis();
        
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
        
        typedef pair<Value*, Value*> LoopBound;
        
        struct LoopInfoStruct{
            vector<Value *>* IDV;                            // induction variable
            vector<Value *>* INC;                            // stride
            vector<LoopIndvBoundAnalysis::LoopBound>* LB;    // bounds
            vector<llvm::CmpInst::Predicate>* PREDICATE;       // upper bound condition >, >=, <, <=
        };
        
        /* The nodes in Loop tree can either be a loop or an array access */
        struct LoopRefTNode {
            int LoopLevel;
            Loop* L;
            LoopInfoStruct* LIS;
            Instruction* AA;
            vector<LoopRefTNode *>* next;
        };
        
        /* Root node for loop reference tree */
        LoopRefTNode* LoopRefTree;
        
        /* Find the Loop Bound */
        LoopBound findLoopBound(Loop *L, Value *var);
        
        /* Find All the Cond/Inc BasicBlocks within the SubLoops */
        vector<BasicBlock *> getSubLoopCondBlock(Loop *L);
        vector<BasicBlock *> getSubLoopIncBlock(Loop *L);
        
        /* Print Func */
        //void dumpLoopInfoStruct();
        
        /* Get the Loop Bound in a string format */
        string getBound(Value * expr);
        
        /* Predicate to string */
        string predicateToString(llvm::CmpInst::Predicate PREDICATE);
        
        void DumpLoopTree(LoopRefTNode* LTroot, std::string prefix);
        LoopInfoStruct* ExtractLoopInfo(Loop *L);
        void LoopTreeConstruction(Loop* LI , LoopRefTNode * root, int level);
        LoopRefTNode* LoopTreeConstructionLoop(LoopRefTNode * root);
        LoopRefTNode* ExtractRefInfo(Instruction* Inst);
        LoopRefTNode* LoopTreeConstructionRef(LoopRefTNode * root, vector<BasicBlock*> BBList);
        
        //void inductionVariableAnalysis(Function &F);
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
        
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}



#endif /* loopAnalysis_hpp */
