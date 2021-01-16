#ifndef uiAccCodeGenOpt_ref_hpp
#define uiAccCodeGenOpt_ref_hpp

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

// #define REFERENCE_GROUP

#define DumpRTMR
// #define DumpRefLease
// #define PSCODEGEN_DEBUG

namespace uiAccCodeGenOpt_ref {
    struct AccLevelUISamplingCodeGenOpt_ref : public FunctionPass {
        static char ID;
        AccLevelUISamplingCodeGenOpt_ref();
        
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
        std::map<Instruction*, std::vector<std::string>> arrayAccessVariable;
        std::vector<Instruction*> outMostIndependentArrayRef;
        
        uint64_t refGlobalNumber = 0;
        std::map<Instruction*, int> refNumber;
        
        std::map<Loop*, int> loopNumber;
        std::map<std::string, int> refToSameArrayCnt;
        std::map<Value*, std::string> indvName;
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> outloops;
        
        std::map<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*, uint64_t> sampleNum;

        void filterArrayAccesses();
        void numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        void initArrayName();
        
        void addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<std::string> indvs);
        void addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        
        void subBlkRTGen(); 
        /* This code used to speedup the reuse search process */
        void rtSearchOptGen();
        void rtHistoGen();
#ifdef DumpRTMR
        void rtDumpGen();
        void rtToMRGen();
        void mrDumpGen();
#elif defined(DumpRefLease)
        void accessRatioCalGen();
        void initHitsCostsGen();
        void getPPUCGen();
        void getMaxPPUCGen();
        void DumpRIGen();
        void RLGen();
#endif
        void headerGen();

        //*********************************************************************
        // Function to handle parallel excution 
        //*********************************************************************
        void initOutLoop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        /* Convert the per-ref rtHist to whole-prog rtHist*/
        void rtMergeGen();

/* Distribute the UI RT based on Gaussian Distribution*/
#if defined(GAUSSIAN_SMOOTHING)
#ifdef REFERENCE_GROUP
        void GroupGaussianDistrGen(string space);
#endif
        void GaussianDistrGen();
/* Distribute the UI RT uniformly */
#elif defined(UNIFORM_SMOOTHING)
#ifdef REFERENCE_GROUP
        void GroupUniformDistrGen(string space);
#endif
        void UniformDistrGen();
#endif // end of SMOOTHING Macro   
        string getBound(Value* bound);
        string getBound_Start(Value* bound);
        string getLoopInc(Value *inc);
        uint64_t getOuterMostLoopIterationSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * node);
        uint64_t computeIterSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        /* 
         * When enableOPT is close
         * - Find loops that contains the sampled referenceID
         * When enableOPT is open
         * - Find loops that contains the reference, the reference could have different ID but should have the same name 
         * The useID here is to filter out those references that will not be iterated anymore because it appears in the loops that had already finished execution 
        */
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> findLoops(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, 
            std::string refName, 
            int useID, 
            std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops,
            bool enableOPT
        );
        
        void searchReuseDifferentLoopsUpdateFuncGen();
        void searchReuseDifferentLoopsCalFuncGen();

        bool searchReuseDifferentLoopsInitGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, 
            bool GenFlag, 
            std::string refName, 
            int useID,
            std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, 
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, 
            string space);

        bool searchReuseDifferentLoopsBodyGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, 
            bool GenFlag, 
            std::string refName, 
            int useID,
            std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, 
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, 
            string space);
        
        void searchReuseSameLoopInitGen(
            std::string refName, 
            int useID, 
            std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, 
            string space);

        void searchReuseSameLoopBodyGen(
            std::string refName, 
            int useID, 
            std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, 
            string space);
        
        bool refRTSearchGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, 
            bool GenFlag, 
            bool isFirstOutLoop, 
            std::string refName, 
            int useID, 
            int sample_number,
            std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, 
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, 
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> sampleIDVs, 
            string space);

        void refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);

        void refRTBodyGen(
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, 
            std::string refName, 
            int useID);
        
        void mainGen();
        
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif
