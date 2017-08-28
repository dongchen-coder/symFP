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

/*
 * Sampling mode 0: No sampling
 * Sampling mode 1: Uniform sampling
 * Sampling mode 2: Random sampling
 */

#define SAMPLING 1

#define UNIFORM_SAMPLING_RATE 0.01

namespace ssCodeGen {
    
    struct StaticSamplingCodeGen : public FunctionPass {
        
        static char ID;
        StaticSamplingCodeGen();
        
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
        
        std::map<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*, int> loopOrder;
        std::map<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*, int> refOrder;
        std::map<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*, int> refCntOfLoop;
        
        std::map<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*, int> refTotalOrder;
        
        std::map<Value*, std::string> indvName;
        
        std::map<std::string, int> refToSameArrayCnt;
        std::map<Instruction*, int> refNumber;
        void numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        void addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<std::string> indvs);
        void addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        void pairRefRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int useID, int reuseID);
        void pairRefRTTopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        void useLoopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops);
        std::string getBound(Value *bound);
        
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* findRef(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int refID);
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> findRefLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int refID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops);
        void reuseLoopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int useID, int reuseID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops);
        
        void checkLocGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int useID, int reuseID);
        
        void checkIntervenRefGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops, int useID, int reuseID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> loops, loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTreeRoot);
        void checkIntervenBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, int useID, int reuseID);
        void checkIntervenTopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        void rtCalFuncBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::string refName, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> useLoops, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> reuseLoops, int useID, int reuseID);
        void rtCalFuncTopGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        void rtHistoGen();
        void rtDumpGen();
        
        void headerGen();
        
        void mainGen();
        
        int initRefTotalOrder(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, int order);
        void initArrayName();
        void initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        int initRefCntOfLoop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        void initRefOrder(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        void initLoopOrder(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        void init(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
        
    };
}

#endif /* ssCodeGen_hpp */
