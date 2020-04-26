#ifndef riCodeGen_ref_hpp
#define riCodeGen_ref_hpp

#include <map>

#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
#include "loopTreeTransform.hpp"
#include "sampleNumAnalysis.hpp"

using namespace llvm;

#define CLS 64
#define DS 8

#define SAMPLING 2

// #define PARALLEL_CXX_THREAD

//#define PROFILE_SEARCH_REUSE

//#define SEARCH_REUSE_SAME_LOOP
//#define SEARCH_REUSE_DIFFERENT_LOOPS

#define DumpRTMR
// #define DumpRefLease

namespace riCodeGen_ref {
    struct AllLevelRICodeGen_ref : public FunctionPass {
        static char ID;
        AllLevelRICodeGen_ref();
        
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
        
        uint64_t refGlobalNumber = 0;
        std::map<Instruction*, int> refNumber;
        
        std::map<Loop*, int> loopNumber;
        std::map<std::string, int> refToSameArrayCnt;
        std::map<Value*, std::string> indvName;
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> outloops;
        
        std::map<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*, uint64_t> sampleNum;
        
        void numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        void initArrayName();
        
        void addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<std::string> indvs);
        void addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        
        void subBlkRTGen(); 
        void rtHistoGen();
#ifdef DumpRTMR
        void rtDumpGen();
        void rtToMRGen();
        void mrDumpGen();
#endif
        void headerGen();

#ifdef PARALLEL
        //*********************************************************************
        // Function to handle parallel excution 
        //*********************************************************************
        /* Convert the per-ref rtHist to whole-prog rtHist*/
        void rtMergeGen();
        /* Distribute the UI RT based on Gaussian Distribution*/
        void GaussianDistrGen();
        /* Distribute the UI RT uniformly */
        void UniformDistrGen();
#endif 
        
        string getBound(Value* bound);
        string getBound_Start(Value* bound);
        
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> findLoops(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree,
            std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> ref
        );

        void LoopIterIncGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* root,
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node,
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops,
            string space
        );
        
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* findFirstRefInLoop(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* outMostLoop
        );
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* findFirstRefAfterLoop(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree,
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* loop
        );
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* findFirstRefAfterRef(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree,
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* Ref
        );
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* findLastRefInLoop(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* loop
        );
        
        bool refRTSearchGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *root,
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree,
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, 
            string space);

        void refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);

        void refRTBodyGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree
        );
        
        void mainGen();
        
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif
