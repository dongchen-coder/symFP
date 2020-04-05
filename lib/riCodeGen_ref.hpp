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
        
        /* 
         * When enableOPT is close
         * - Find loops that contains the sampled referenceID
         * When enableOPT is open
         * - Find loops that contains the reference, the reference could have different ID but should have the same name 
         * The useID here is to filter out those references that will not be iterated anymore because it appears in the loops that had already finished execution 
        */
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> findLoops(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, 
            //std::string refName,
            //int useID,
            std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops,
            bool enableOPT
        );

        void LoopIterIncGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node,
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops,
            string space
        );
        
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* findFirstRef(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* outMostLoop
        );
        
        bool refRTSearchGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree,
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, 
            string space);

        void refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);

        void refRTBodyGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree
            //std::string refName,
            //int useID
        );
        
        void mainGen();
        
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif
