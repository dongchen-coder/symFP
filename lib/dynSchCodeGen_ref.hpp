#ifndef dynSchCodeGen_ref_hpp
#define dynSchCodeGen_ref_hpp

#include <map>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
#include "loopTreeTransform.hpp"
#include "IVDependenceAnalysis.hpp"
#include "AccessGraphAnalysis.hpp"
#include "plumCodeGenUtil.hpp"
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

typedef std::pair<DepNode *, int> DepParentInfo;

namespace dynSchCodeGen_ref {
    struct DynamicSchedulingSamplerCodeGen : public FunctionPass {
        static char ID;
        DynamicSchedulingSamplerCodeGen();
        
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
		unordered_set<std::string> ArrayNameSet;
        
        uint64_t refGlobalNumber = 0;
        std::map<Instruction*, int> refNumber;
        
        std::map<Loop*, int> loopNumber;
        std::map<std::string, int> refToSameArrayCnt;
        std::map<Value*, std::string> indvName;
        std::vector<LoopRefTNode*> outloops;
        
        std::map<LoopRefTNode*, uint64_t> sampleNum;
		
		double SamplingRate = 0.0;
        
        void numberRefToSameArray(LoopRefTNode *LoopRefTree);
        void numberLoops(LoopRefTNode *LoopRefTree);
        void initIndvName(LoopRefTNode* LoopRefTree);
        void initArrayName();
        
        void addrCalFuncGen(LoopRefTNode* LoopRefTree, std::vector<std::string> indvs);
        void addrCalFuncGenTop(LoopRefTNode *LoopRefTree);
        
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
        
		std::string getBound(Value* bound);
		std::string getBound_Start(Value* bound);
		DepParentInfo getParent(Value * v);
        
        std::vector<LoopRefTNode*> findLoops(
            LoopRefTNode *LoopRefTree,
            std::vector<LoopRefTNode*> ref
        );

		void LoopIterUpdateGen(AccGraph *G,
							   AccGraphEdge *E,
							   LoopRefTNode* IDomLoop,
							   std::string space,
							   bool EnableIterationInc,
							   bool EnableIterationUpdate);
		
		void reuseUpdateGen(LoopRefTNode * access,
							int currentAccessLoopNestCnt,
							std::string accessName,
							std::string space);
        
        bool refRTSearchGen(LoopRefTNode *root,
							LoopRefTNode *LoopRefTree,
							std::vector<LoopRefTNode*> currentLoops,
							std::string space);
		
		bool perArrayRTSearchGen(std::string array,
								 LoopRefTNode *root,
								 LoopRefTNode *LoopRefTree,
								 std::vector<LoopRefTNode*> currentLoops,
								 std::string space);
		
		void perArrayRTBodyGen(std::string array, LoopRefTNode *LoopRefTree);

        void refRTGen(LoopRefTNode *LoopRefTree);

        void refRTBodyGen(
            LoopRefTNode *LoopRefTree
        );
        
        void mainGen();
        
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif /* dynSchCodeGen_hpp */
