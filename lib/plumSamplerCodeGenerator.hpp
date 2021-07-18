//
//  plumSamplerCodeGenerator.hpp
//  LLVMSPS
//
//  Created by noya-fangzhou on 7/13/21.
//

#ifndef plumSamplerCodeGenerator_hpp
#define plumSamplerCodeGenerator_hpp

#include <string>
#include <vector>
#include <map>
#include "loopAnalysis.hpp"
#include "IVDependenceAnalysis.hpp"
#include "plumCodeGenUtil.hpp"
#include "AccessGraphAnalysis.hpp"

using namespace std;

enum SchedulingType {
	STATIC,
	DYNAMIC
};

typedef std::pair<DepNode *, int> DepParentInfo;
class SamplerCodeGenerator {
public:
	SamplerCodeGenerator() {}
	SamplerCodeGenerator(loopAnalysis::LoopIndvBoundAnalysis * LA,
						 ivdepAnalysis::IVDependenceAnalysis * IVA,
						 AccGraphAnalysis::AccessGraphAnalysis * GA,
						 std::map<Instruction*, int> refNumber) : LA(LA), IVA(IVA), GA(GA), refNumber(refNumber) {}
	void setSamplingRate(double rate) {
		this->SamplingRate = rate;
	}
	void setSchedulingType(SchedulingType newType) {
		this->type = newType;
	}
	void rtHistoGen();
	/* generate rtToMR function */
	void rtToMRGen();
	/* generate rtDump function */
	void rtDumpGen();
	/* generate mrDump function */
	void mrDumpGen();
	
	DepParentInfo getParent(Value * V);
	
	void RTBodyGen(LoopRefTNode *LoopRefTree);
	void perArrayRTBodyGen(std::string array, LoopRefTNode *LoopRefTree);
	
	void refRTSearchGen(LoopRefTNode *root,
						LoopRefTNode *LoopRefTree,
						std::vector<LoopRefTNode*> & outloops,
						std::string space);
	
	void perArrayRTSearchGen(std::string array,
							 LoopRefTNode *root,
							 LoopRefTNode *LoopRefTree,
							 std::vector<LoopRefTNode*> & outloops,
							 std::string space);

private:
	SchedulingType type = STATIC; // default is static scheduling
	double SamplingRate;
	loopAnalysis::LoopIndvBoundAnalysis * LA;
	ivdepAnalysis::IVDependenceAnalysis * IVA;
	AccGraphAnalysis::AccessGraphAnalysis * GA;
	std::map<Instruction*, int> refNumber;
	
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
	
	void StaticRefRTSearchGen(LoopRefTNode *root,
						LoopRefTNode *LoopRefTree,
						std::vector<LoopRefTNode*> & outloops,
						std::string space);
	
	void StaticPerArrayRTSearchGen(std::string array,
							 LoopRefTNode *root,
							 LoopRefTNode *LoopRefTree,
							 std::vector<LoopRefTNode*> & outloops,
							 std::string space);
	
	void DynamicRefRTSearchGen(LoopRefTNode *root,
						LoopRefTNode *LoopRefTree,
						std::vector<LoopRefTNode*> & outloops,
						std::string space);
	
	void DynamicPerArrayRTSearchGen(std::string array,
							 LoopRefTNode *root,
							 LoopRefTNode *LoopRefTree,
							 std::vector<LoopRefTNode*> & outloops,
							 std::string space);
};

#endif /* plumSamplerCodeGenerator_hpp */
