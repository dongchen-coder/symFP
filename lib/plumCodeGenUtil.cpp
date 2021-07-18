//
//  plumCodeGenUtil.cpp
//  LLVMSPS
//
//  Created by noya-fangzhou on 5/11/21.
//

#include "plumCodeGenUtil.hpp"
#include <algorithm>

#define DEBUG_TYPE	"plum-util"

using namespace std;

AccessQueryResult::AccessQueryResult() {}
AccessQueryResult::AccessQueryResult(LoopRefTNode * anode, vector<LoopRefTNode *> loops)
{
	this->AccessNode = anode;
	this->LoopsToAccess = loops;
}

AccessQueryEngine::AccessQueryEngine() {}
AccessQueryEngine::AccessQueryEngine(LoopRefTNode * root)
{
	assert(!root->AA);
	this->root = root;
	buildParentLoopRefNodeTree(root);
	this->accessIter = this->AccessNodeSet.begin();
}
void AccessQueryEngine::buildParentLoopRefNodeTree(LoopRefTNode * root)
{
	// iterate the LoopRefTree using DFS
	if (root->AA) { // access node
		this->AccessNodeSet.emplace_back(root);
	} else if (root->L) { // loop node
		for (auto subNode : *(root->next)) {
			this->ParentLoopRefTNodeMap.insert(make_pair(subNode, root));
			buildParentLoopRefNodeTree(subNode);
		}
		LoopRefTNode * lastAccessNode = *(root->next->rbegin());
		if (lastAccessNode->AA) {
			// we store those access nodes that are also the last child node
			// of root.
			//
			// Only store the access node (A) in case 2
			// case 1: L -> ... -> A -> L
			// case 2: L -> ... -> L -> A
			this->LastAccessNodeInLoop.insert(lastAccessNode);
		}
	}
	return;
}
void AccessQueryEngine::getPath(LoopRefTNode * node, LoopRefTNode * terminatePoint, vector<LoopRefTNode *> &path)
{
	LoopRefTNode * tmpNode = node;
	while (this->ParentLoopRefTNodeMap.find(tmpNode) != this->ParentLoopRefTNodeMap.end()) {
		if (this->ParentLoopRefTNodeMap[tmpNode] != terminatePoint) {
			path.emplace_back(this->ParentLoopRefTNodeMap[tmpNode]);
			tmpNode = this->ParentLoopRefTNodeMap[tmpNode];
		} else {
			break;
		}
	}
}
vector<LoopRefTNode *> AccessQueryEngine::getAccessNodeInOrder()
{
	return this->AccessNodeSet;
}
bool AccessQueryEngine::hasNext()
{
	return this->accessIter != this->AccessNodeSet.end();
}
bool AccessQueryEngine::isFirstAccess(LoopRefTNode * node)
{
	return find(this->AccessNodeSet.begin(), this->AccessNodeSet.end(), node) == this->AccessNodeSet.begin();
}

bool AccessQueryEngine::isLastAccessNodeInALoop(LoopRefTNode *node)
{
	return this->LastAccessNodeInLoop.find(node) != this->LastAccessNodeInLoop.end();
}

bool AccessQueryEngine::areTwoAccessInSameLevel(LoopRefTNode * accA, LoopRefTNode * accB)
{
	LoopRefTNode * parentA = NULL, * parentB = NULL;
	if (this->ParentLoopRefTNodeMap.find(accA) != this->ParentLoopRefTNodeMap.end()) {
		parentA = this->ParentLoopRefTNodeMap[accA];
	}
	if (this->ParentLoopRefTNodeMap.find(accB) != this->ParentLoopRefTNodeMap.end()) {
		parentB = this->ParentLoopRefTNodeMap[accB];
	}
	if ((!parentA && parentB) || (parentA && !parentB))
		return false;
	return parentA->LoopLevel == parentB->LoopLevel || (!parentA && !parentB);
}

bool AccessQueryEngine::areTwoAccessInSameLoop(LoopRefTNode * accA, LoopRefTNode * accB)
{
	LoopRefTNode * parentA = NULL, * parentB = NULL;
	if (this->ParentLoopRefTNodeMap.find(accA) != this->ParentLoopRefTNodeMap.end()) {
		parentA = this->ParentLoopRefTNodeMap[accA];
	} else {
		return false;
	}
	if (this->ParentLoopRefTNodeMap.find(accB) != this->ParentLoopRefTNodeMap.end()) {
		parentB = this->ParentLoopRefTNodeMap[accB];
	} else {
		return false;
	}
	return parentA == parentB;
}
AccessQueryResult * AccessQueryEngine::getNext()
{
	if (!hasNext()) {
		this->accessIter = this->AccessNodeSet.begin();
		return NULL; // if no acccesses any more, return NULL
	}
	vector<LoopRefTNode *> loops; // store list to the root in reverse order, from inner to outer
	LoopRefTNode * tmpNode = *(this->accessIter);
	LoopRefTNode * nextAccess = *(this->accessIter);
	this->getPath(tmpNode, NULL, loops);
	this->accessIter++; // proceed to the next
	LLVM_DEBUG(dbgs() << " has next " <<
			   *(nextAccess->AA) << "\n");
	return new AccessQueryResult(nextAccess, loops);
	
}
LoopRefTNode * AccessQueryEngine::getFirstAccessInLoop(LoopRefTNode *node)
{
	if (node->L) { // this is a loop node, we dig one step deeper
		return getFirstAccessInLoop(node->next->front());
	} else if (node->AA) {
		return node;
	}
	return NULL;
}

LoopRefTNode * AccessQueryEngine::getLastAccessInLoop(LoopRefTNode *node)
{
	if (node->L) { // this is a loop node
		// get the last node in node->next, we call it last_child_node
		// if last_child_node is a access node, then this is the last access node
		// if last_child_node is a loop node, we recursively find its last access node
		LoopRefTNode * lastChildNode = *(node->next->rbegin());
		if (lastChildNode->L) {
			return getLastAccessInLoop(lastChildNode);
		} else if (lastChildNode->AA) {
			return lastChildNode;
		}
//		else {
//			 assert(true, "ERROR! This node is neither a access nor a loop node!");
//		}
	} else if (node->AA) {
		return node;
	}
	return NULL;
}

void AccessQueryEngine::getImmediateCommonLoopDominator(LoopRefTNode * nodeA, LoopRefTNode * nodeB, bool isLastAccess, vector<LoopRefTNode *> & dominators)
{
	vector<LoopRefTNode *> pathToA, pathToB;
	this->getPath(nodeA, NULL, pathToA);
	this->getPath(nodeB, NULL, pathToB);
	vector<LoopRefTNode *>::iterator pathAIter = pathToA.begin();
	vector<LoopRefTNode *>::iterator pathBIter = pathToB.begin();
	int degreeDiff = pathToB.size() - pathToA.size();
	while (degreeDiff != 0) {
		if (degreeDiff > 0) { // B is deeper than A
			pathBIter++;
			degreeDiff--;
		} else if (degreeDiff < 0) { // A is deeper than B
			pathAIter++;
			degreeDiff++;
		}
	}
	// now pathAIter and pathBIter at the same loop level
	while (pathAIter != pathToA.end() && pathBIter != pathToB.end()) {
		if ((*pathAIter) == (*pathBIter)) { // these two loops are the same
			dominators.emplace_back(*pathAIter);
			break;
		}
		pathAIter++;
		pathBIter++;
	}
	// all access nodes are dominated by the outermost loop
	// dominators.size() should >= 1
	assert(dominators.size() == 1);
}
LoopRefTNode * AccessQueryEngine::getImmediateLoopDominator(LoopRefTNode * node)
{
	if (this->ParentLoopRefTNodeMap.find(node) != this->ParentLoopRefTNodeMap.end()) {
		return this->ParentLoopRefTNodeMap[node];
	}
	return NULL;
}
string AccessQueryEngine::toString()
{
	return "";
}

/* AccessGraph Analysis */
AccGraphEdge::AccGraphEdge(LoopRefTNode * source, LoopRefTNode * sink, LoopRefTNode * loop)
{
	this->source = source;
	this->sink = sink;
	this->carryLoop = loop;
}

bool AccGraphEdge::containsNode(LoopRefTNode * node)
{
	assert(node->AA);
	return this->source == node || this->sink == node;
}

bool AccGraphEdge::hasSource(LoopRefTNode * node)
{
	assert(node->AA);
	return this->source == node;
}

bool AccGraphEdge::hasSink(LoopRefTNode * node)
{
	assert(node->AA);
	return this->sink == node;
}

LoopRefTNode * AccGraphEdge::getSourceNode()
{
	return this->source;
}

LoopRefTNode * AccGraphEdge::getSinkNode()
{
	return this->sink;
}

LoopRefTNode * AccGraphEdge::getCarryLoop()
{
	return this->carryLoop;
}

bool AccGraphEdge::operator==(const AccGraphEdge & rhs)
{
	return ((this->source == rhs.source && this->sink == rhs.sink) || (this->source == rhs.sink && this->sink == rhs.source)) && this->carryLoop == rhs.carryLoop;
}

bool AccGraphEdge::operator<( const AccGraphEdge & rhs) const
{
	if (this->carryLoop && !rhs.carryLoop) // edge with carryLoop is smaller
		return true;
	else if (this->carryLoop && rhs.carryLoop) // two edges have carryLoop, lower caried loop level is smaller
		return this->carryLoop->LoopLevel < rhs.carryLoop->LoopLevel;
	else // both edge with no carryLoop, any order is okay. Will not have because An AccessGraphNode can have one and only one normal edge
		return true;
	
}

AccGraph::AccGraph() {}
AccGraph::AccGraph(LoopRefTNode * root)
{
	this->QueryEngine = AccessQueryEngine(root);
	buildAccGraph();
}
void AccGraph::buildAccGraph()
{
	AccessQueryResult * current = NULL, * next = NULL;
	current = this->QueryEngine.getNext();
	while (true) {
		if (!this->root) {
			// assign it to the root if it is NULL
			this->root = current->AccessNode;
		}
		this->GraphNodeSet.insert(current->AccessNode);
		next = this->QueryEngine.getNext();
		// current access node is also the last node of a loop
		// in this case there is a possibility that a backward edge exists
		// we will examine it with extra steps
		if (QueryEngine.isLastAccessNodeInALoop(current->AccessNode)) {
			// iterate all previous loop nodes, from inner to outer
			// for each loop, we found its first access node.
			// Then there exists a backedge that connect current with this first access
			// node.
			//
			// Note:
			// 1. if next is not NULL && path to current only has 1 loop (outermost)
			//		we do not add any backedges
			// 2. if next && path to curent contains multiple loops, we build backedge
			//		for those inner loops (except the last one)
			// 3. if !next, we build backedge with the outermost loop
			if (!next || current->LoopsToAccess.size() > 1) { // case 2, 3
				vector<LoopRefTNode *>::iterator it = current->LoopsToAccess.begin();
				LoopRefTNode * firstAA = NULL;
				while (it != current->LoopsToAccess.end()) {
					if (next && it == current->LoopsToAccess.end()-1) { // case 2
						break;
					}
					if (QueryEngine.getLastAccessInLoop(*it) == current->AccessNode) {
						firstAA = QueryEngine.getFirstAccessInLoop(*it);
						assert(firstAA && this->GraphNodeSet.find(firstAA) != this->GraphNodeSet.end());
						AccGraphEdge * backEdge = new AccGraphEdge(current->AccessNode, firstAA, *it);
						this->GraphEdgeSet.emplace_back(backEdge);
					}
					it++;
				}
			}
		}
		if (next) {
			AccGraphEdge * edge = new AccGraphEdge(current->AccessNode, next->AccessNode, NULL);
			this->GraphEdgeSet.emplace_back(edge);
			current = next;
			continue;
		}
		// only if current is the last access node will reach here
		break;
	}
}
LoopRefTNode * AccGraph::getRoot()
{
	return this->root;
}
vector<LoopRefTNode *> AccGraph::getAllGraphNodes()
{
	return this->QueryEngine.getAccessNodeInOrder();
}
vector<AccGraphEdge *> AccGraph::getAllGraphEdges()
{
	return this->GraphEdgeSet;
}
vector<AccGraphEdge *> AccGraph::getEdgesWithSource(LoopRefTNode * source)
{
	vector<AccGraphEdge *> edges;
	for (auto Edge : this->GraphEdgeSet) {
		if (Edge->hasSource(source))
			edges.emplace_back(Edge);
	}
	return edges;
//	// if multiple edges are found
//	// sort them based on the following order:
//	// 1. those edges have carried loop, we sort them based on the carry-loop's level
//	// lower loop level (outer loop) will be iterated first
//	// 2. edge with no carry loop
//	sortEdges(edges);
//	return edges;
}

void AccGraph::sortEdges(std::vector<AccGraphEdge *> & edges)
{
	if (edges.size() == 1)
		return;
	std::sort(edges.begin(), edges.end());
}

vector<AccGraphEdge *> AccGraph::getEdgesWithSink(LoopRefTNode * sink)
{
	vector<AccGraphEdge *> edges;
	for (auto Edge : this->GraphEdgeSet) {
		if (Edge->hasSink(sink))
			edges.emplace_back(Edge);
	}
	return edges;
}

vector<LoopRefTNode *> AccGraph::getPathToNode(LoopRefTNode * from, LoopRefTNode * to)
{
	vector<LoopRefTNode *> path;
	this->QueryEngine.getPath(from, to, path);
	return path;
}

/* IDV Dependence Analysis */
DepNode::DepNode(int level, LoopInfoStruct * v)
{
   this->level = level;
   this->IV = v;
}

bool DepNode::isParent(DepNode * node)
{
   return this->Parents.contains(node);
}

bool DepNode::isDependence(DepNode * node)
{
   return this->Dependences.contains(node);
}

void DepNode::insertParent(DepNode * node)
{
   this->Parents.insert(node);
}

void DepNode::insertDependence(DepNode * node)
{
   this->Dependences.insert(node);
}

bool DepNode::hasDependence()
{
   return this->Dependences.empty();
}

string DepNode::toString()
{
   return "";
}
