//
//  plumCodeGenUtil.hpp
//  LLVMSPS
//
//  Created by noya-fangzhou on 5/11/21.
//

#ifndef plumCodeGenUtil_hpp
#define plumCodeGenUtil_hpp

#include <stdio.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <vector>
#include <map>
#include "llvm/ADT/SmallSet.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/DenseMap.h"
#include "loopAnalysis.hpp"

using namespace llvm;

typedef loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode LoopRefTNode;
typedef loopAnalysis::LoopIndvBoundAnalysis::LoopInfoStruct LoopInfoStruct;

/* This engine return the access and the loops it traversed */
class AccessQueryResult {
public:
	LoopRefTNode * AccessNode;
	std::vector<LoopRefTNode *> LoopsToAccess;
	AccessQueryResult();
	AccessQueryResult(LoopRefTNode * anode, vector<LoopRefTNode *> loops);
};
class AccessQueryEngine {
private:
	LoopRefTNode * root;
	std::unordered_map<LoopRefTNode*, LoopRefTNode*> ParentLoopRefTNodeMap;
	std::unordered_set<LoopRefTNode*> LastAccessNodeInLoop; /* store the set of access nodes that are the last access in a loop */
	std::vector<LoopRefTNode *> AccessNodeSet; /* store the set of access node based on the traverse order */
	void buildParentLoopRefNodeTree(LoopRefTNode * root);
public:
	std::vector<LoopRefTNode *>::iterator accessIter;
	AccessQueryEngine();
	AccessQueryEngine(LoopRefTNode * root);
	bool hasNext();
	bool isFirstAccess(LoopRefTNode * node);
	bool isLastAccessNodeInALoop(LoopRefTNode * node);
	bool areTwoAccessInSameLoop(LoopRefTNode * accA, LoopRefTNode * accB);
	bool areTwoAccessInSameLevel(LoopRefTNode * accA, LoopRefTNode * accB);
	vector<LoopRefTNode *> getAccessNodeInOrder();
	void getPath(LoopRefTNode * node, LoopRefTNode * terminatePoint, vector<LoopRefTNode *> &path);
	LoopRefTNode * getFirstAccessInLoop(LoopRefTNode * node);
	LoopRefTNode * getLastAccessInLoop(LoopRefTNode *node);
	// get the common loop dominator of two accesses
	// the common loop dominator does contain the root outermost loop
	// SPECIAL CASE: if is the last access node of the entire outermost loops,
	// then we need to consider the outermost loop
	void getImmediateCommonLoopDominator(LoopRefTNode * nodeA, LoopRefTNode * nodeB, bool isLastAccess, std::vector<LoopRefTNode *> & dominators);
	AccessQueryResult * getNext();
	LoopRefTNode * getImmediateLoopDominator(LoopRefTNode * node);
	string toString();
};

/* AccessGraph Analysis */
class AccGraphEdge
{
private:
	LoopRefTNode * source;
	LoopRefTNode * sink;
	LoopRefTNode * carryLoop; // NULL means forward edge else backward edge
public:
	AccGraphEdge();
	AccGraphEdge(LoopRefTNode * source, LoopRefTNode * sink, LoopRefTNode * loop);
	std::string toString();
	LoopRefTNode * getCarryLoop();
	LoopRefTNode * getSourceNode();
	LoopRefTNode * getSinkNode();
	bool containsNode(LoopRefTNode * node);
	bool hasSource(LoopRefTNode * node);
	bool hasSink(LoopRefTNode * node);
	bool operator==(const AccGraphEdge & rhs );
	bool operator<( const AccGraphEdge & rhs) const;
};
class AccGraph
{
private:
	std::set<LoopRefTNode *> GraphNodeSet;
	std::vector<AccGraphEdge *> GraphEdgeSet;
	LoopRefTNode * root;
	void buildAccGraph();
	void buildAccGraph(std::string arrayName); // build the access graph for all access to specific array
public:
	AccessQueryEngine QueryEngine;
	AccGraph();
	AccGraph(LoopRefTNode * root);
	AccGraph(LoopRefTNode * root, std::string arrayName);
	LoopRefTNode * getRoot();
	std::vector<LoopRefTNode *> getAllGraphNodes();
	std::vector<LoopRefTNode *> getPathToNode(LoopRefTNode * from, LoopRefTNode * to);
	std::vector<AccGraphEdge *> getAllGraphEdges();
	std::vector<AccGraphEdge *> getEdgesWithSource(LoopRefTNode * source);
	std::vector<AccGraphEdge *> getEdgesWithSink(LoopRefTNode * sink);
private:
	void sortEdges(std::vector<AccGraphEdge *> & edges);
};



/* IDV Dependence Analysis */
class DepNode
{
public:
	int level;
	LoopInfoStruct * IV;
	SmallSet<DepNode *, 5> Dependences;
	SmallSet<DepNode *, 5> Parents;

	DepNode(int level, LoopInfoStruct * v);
	~DepNode();
	bool isParent(DepNode * node);
	bool isDependence(DepNode * node);
	bool hasDependence();
	void insertParent(DepNode * node);
	void insertDependence(DepNode * node);
	string toString();
};



#endif /* plumCodeGenUtil_hpp */
