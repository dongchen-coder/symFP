//
//  plumCodeGen_ref.hpp
//  LLVMSPS
//
//  Created by noya-fangzhou on 7/13/21.
//

#ifndef plumCodeGen_ref_hpp
#define plumCodeGen_ref_hpp

#include <map>
#include <unordered_set>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
#include "loopTreeTransform.hpp"
#include "IVDependenceAnalysis.hpp"
#include "AccessGraphAnalysis.hpp"
#include "plumSamplerCodeGenerator.hpp"
#include "plumCodeGenUtil.hpp"
#include "sampleNumAnalysis.hpp"

using namespace llvm;

#define CLS 64
#define DS 8

#define SAMPLING 2
// #define PER_ARRAY_ANALYSIS

// #define PARALLEL_CXX_THREAD

//#define PROFILEc_SEARCH_REUSE

// #define DumpRefLease

namespace plumCodeGen_ref {
	struct PLUMSamplerCodeGen : public FunctionPass {
		static char ID;
		PLUMSamplerCodeGen();

		std::map<Instruction*, std::string> arrayExpression;
		unordered_set<std::string> ArrayNameSet;
		
		uint64_t refGlobalNumber = 0;
		std::map<Instruction*, int> refNumber;
		std::map<std::string, int> refToSameArrayCnt;
		std::map<Value*, std::string> indvName;
		
		std::vector<LoopRefTNode*> outloops;
		
		SamplerCodeGenerator * CodeGen;
		
		loopAnalysis::LoopIndvBoundAnalysis * LA;
		ivdepAnalysis::IVDependenceAnalysis * IVA;
		AccGraphAnalysis::AccessGraphAnalysis * GA;
		
		
		double SamplingRate = 0.0;
		
		void initIndvName(LoopRefTNode* LoopRefTree);
		void numberRefToSameArray(LoopRefTNode *LoopRefTree);
		
		void addrCalFuncGen(LoopRefTNode* LoopRefTree, std::vector<std::string> indvs);
		void addrCalFuncGenTop(LoopRefTNode *LoopRefTree);

		void headerGen();
		
		void refRTGen(LoopRefTNode *LoopRefTree);
		
		void perArrayRTBodyGen(std::string array, LoopRefTNode *LoopRefTree);
		void refRTBodyGen(LoopRefTNode *LoopRefTree);
		
		void mainGen();
		
		bool runOnFunction(Function &F) override;
		void getAnalysisUsage(AnalysisUsage &AU) const override;
	};
}

#endif /* plumCodeGen_ref_hpp */
