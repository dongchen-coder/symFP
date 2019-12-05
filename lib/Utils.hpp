#ifndef Utils_hpp
#define Utils_hpp

#include <vector>

using namespace std;

/* The base class of All LoopNode */
struct TreeNodeBase
{
public:
	vector<TreeNodeBase *>* next;
	TreeNodeBase();
	~TreeNodeBase();
	
};

#endif