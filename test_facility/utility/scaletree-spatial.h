#include <iostream>
#include <fstream>
#include "math.h"
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"

#include "histo.H"
//#include "fullHisto.h"

#define CandidateSize 1024*1024
#define BUFFERSIZE 4*1024*1024
#define NumCounters 2000
#define LogLinearSize 8
#define HashSize 50000000
#define LOG_BLOCK_SIZE 5
#define REDUCED_THRESHOLD 3

typedef struct tree_node Tree;
struct tree_node {
	Tree * left, * right;
	unsigned long long item;    /* last access time of this node */
	unsigned nodeWt;  /* weight of the node */
	unsigned weight;  /* weight of the entire subtree, including this node */
	unsigned maxSize;    /* maximal size of the node */
	Tree *prev, *next; /* nodes that are next to this node in the sorted order */
};

typedef struct hashEntry {
	unsigned addr;
	unsigned long long cycle;  /* the time of the last access to addr */
	struct hashEntry *hshNxt;
} HashEntry;

using namespace histo;
static const  uint32_t              SUBLOG_BITS = 8;
static const  uint32_t              HIST_BUCKETS = (65-SUBLOG_BITS)*(1<<SUBLOG_BITS);
static histogram<HIST_BUCKETS, sublog_value_to_index<HIST_BUCKETS, SUBLOG_BITS>, sublog_index_to_value<HIST_BUCKETS, SUBLOG_BITS> >  
                             reuse_dist_hist;


/*
inline void RecordDistance(uint64_t dis) {
    reuse_dist_hist.put_value(dis);
//	record(dis);
}
*/

#define RecordDistance(dis)\
{\
	int li;\
	if ((dis)==0)\
	{\
		counters[0]++;\
	}\
	else\
	{\
		li = (dis) >> LogLinearSize;\
		if (li > 0)\
		{\
			if(LogLinearSize+li<NumCounters)\
			{\
				counters[LogLinearSize+li]++;\
			}\
			else\
			{\
				counters[NumCounters-1]++;\
			}\
		}\
		else\
		{\
			int t=0, td=(dis);\
			assert(li==0);\
			while (td>0)\
			{\
				t++;\
				td = td>>1;\
			}\
			counters[t]++;\
		}\
	}\
}


#define MemAddr void *

class ScaleTree {
 public:
    //ScaleTree(int alloc_buf);
    Tree *trace;
    int numInstr;
    int logblocksize;
    int buffersize;
    int logTrace;
    int instr;
    float errorRate;
    unsigned long *buffer;
	unsigned *rd2;
    void processBuffer();
    int PrintResults();
    unsigned long long returnItem(Tree *);
    int bufferindex;
    
    unsigned numData;
    unsigned long long curCycle;
    unsigned power;
    long buffersprocessed;
    HashEntry **hash;
    HashEntry * candidate;
    int hashUsed; 
    unsigned hashallocated;
    unsigned long long counters[NumCounters];
    unsigned sizes[NumCounters];
    float slq[NumCounters];
    int linearSize;    /* size of the linear scale section */
    unsigned int sizeTrace;
    Tree * freeNodeList;
    int Init();
    int DataAccess(unsigned addr);
    void PrintSize();
    int CounterInitialize();
    int HashInitialize();
    unsigned long long HashSearchUpdate(unsigned addr, unsigned long long cyc);
    Tree  * MallocNode();
    void FreeNode(Tree *node);
    Tree * ScaleTreeSplay(unsigned long long i, Tree * t, unsigned *dis);
    Tree * ScaleTreeInsertAtFront(unsigned long long blockEnd, Tree * t, Tree *newNode);
    Tree * QueryScaleTree(unsigned long long oldCyc, unsigned long long newCyc, Tree *t, unsigned *dis);
    Tree *CompactScaleTree(Tree *root);
    void allocateCandidate();
};


int ScaleTree::DataAccess(unsigned addr) {
    unsigned long long lastAccCyc;
    unsigned dis = 0;
    lastAccCyc = HashSearchUpdate(addr, curCycle);
    trace = QueryScaleTree(lastAccCyc, curCycle, trace, &dis);
    if (lastAccCyc == curCycle) {
        /* a new element */
        numData ++;
        if (numData-((numData>>10)<<10) == 0) {
            power=1;
            unsigned tmp = numData;
            while (tmp > 0) {
                tmp = tmp >> 6;
                power ++;
            }
        }
    } else {
        RecordDistance(dis);
    }
    
    //check after every 65K accesses
    if(((sizeTrace & 0xffff) == 0) && (sizeTrace > (power<<13)))
        //compress at size 4*log_{1/errorRare}^M + 4
        //see DingZhong:PLDI03, pp 247
        if(sizeTrace > 4*log((double)numData)/log((double)1.0/(1-errorRate))+4)
            trace = CompactScaleTree(trace);
    
    curCycle ++;
    return dis;
}

void ScaleTree::PrintSize() {
    printf("hash allocated %u,  hash table size %u; scaletree size %u.\n",hashallocated, numData, sizeTrace);
}

int ScaleTree::CounterInitialize() {
    int i;
    int base = 2, logSize = LogLinearSize, cnt=1, log=2;
    sizes[0] = 0;
    while (cnt<=logSize) {
        sizes[cnt] = log;
        log *= base;
        cnt ++;
    }
    linearSize = sizes[logSize];
    
    for (i=logSize+1; i<NumCounters; i++)
        sizes[i] = sizes[i-1] + linearSize;
    
    for (i=0; i<NumCounters; i++) {
        counters[i] = 0;
        slq[i] = 0;
    }
    
    return 0;
}

int ScaleTree::PrintResults() {
    
	int i, lastNonZero=0;
    unsigned long long totAcc=0;
    trace = CompactScaleTree(trace);
    
    std::cout << "/* Total number of data blocks is " << numData << " */" << std::endl;

//  reuse_dist_hist.print();
//	printRDHist();
	
	
	totAcc = numData;
    for (i=0; i<NumCounters; i++) {
        if (counters[i]==0) continue;
        if (i>lastNonZero) lastNonZero = i;
        totAcc += counters[i];
    }

    std::cout << "/* Total number of accesses is " << totAcc << " */" << std::endl;
    std::cout << "/* Bin 0 Distance 0 */\tReuses: "<< counters[0] << std::endl;
    std::cout << "/* Bin 1 Distance 1 */\tReuses: "<< counters[1] << std::endl;

   
    for (i=2; i<=lastNonZero; i++) {
        std::cout << "/* Bin " << i << " Distance " << sizes[i-1] << " to " << sizes[i]-1 <<" */\tReuses: " << counters[i];
        if ((i >= REDUCED_THRESHOLD) && (counters[i] > 0))
            std::cout << "\tSLQ: " << slq[i];
        std::cout << std::endl;
    }
    std::cout << "/* End tree size is " << sizeTrace << " */" << std::endl;
    return lastNonZero;
//	return 0;
}


int ScaleTree::Init() {
    if (sizeof(char *) != sizeof(unsigned long)) {
        printf("the length of a pointer is different with the size of an unsinged long. You need to fix this problem\n");
        exit(0);
    }
    HashInitialize();
    CounterInitialize();
    return 0;
}

int ScaleTree::HashInitialize() {
    int i;
    allocateCandidate();
    hash=(HashEntry**)malloc(sizeof(HashEntry*)*HashSize);
    for (i=0; i<HashSize; i++) hash[i]=NULL;
    return 0;
}

void ScaleTree::allocateCandidate() {
    candidate = (HashEntry*)malloc(sizeof(HashEntry)*CandidateSize);
    if (candidate == NULL) {
        trace = CompactScaleTree(trace);
        printf("Tree is compacted because of insufficient memory for hash\n");
    }
    hashUsed=0;
}

/* addr is accessed at cyc, Cycle cyc means never accessed before
 * insert addr at the head if it is not found
 */
unsigned long long ScaleTree::HashSearchUpdate(unsigned addr, unsigned long long cyc) {
    unsigned hshKey;
    HashEntry *entry;
    unsigned long long oldCyc;
    hshKey = addr % HashSize;
    entry = hash[hshKey];
    
    while (entry != NULL && entry->addr != addr) {
        entry=entry->hshNxt;
    }
    if (entry != NULL) {
        oldCyc = entry->cycle;
        entry->cycle = cyc;
        return oldCyc;
    } else {
        entry = &candidate[hashUsed++];
        if (hashUsed==CandidateSize)
            allocateCandidate();
        entry->addr = addr;
        entry->cycle = cyc;
        entry->hshNxt = hash[hshKey];
        hash[hshKey] = entry;
        return cyc;  /* means that a previous record is not found */
    }
}


Tree *ScaleTree::MallocNode() {
    if (!freeNodeList) {
        Tree * l = (Tree *) malloc(sizeof(Tree)*CandidateSize);
        if (l==NULL) {
            trace = CompactScaleTree(trace);
            printf("Tree is compacted because of insufficient memory for scaletree\n");
        } else {
            int i;
            for (i=0;i<CandidateSize-1;i++)
                l[i].next = &(l[i+1]);
            l[CandidateSize-1].next=freeNodeList;
            freeNodeList = l;
        }
    }
    Tree * t=freeNodeList;
    freeNodeList=freeNodeList->next;
    return t;
}

void ScaleTree::FreeNode(Tree * node) {
    node->next=freeNodeList;
    freeNodeList = node;
}


Tree *ScaleTree::ScaleTreeSplay(unsigned long long i, Tree * t, unsigned *dis) {
    Tree N, *l, *r, *y;
    unsigned left = 0, right = 0;
    if (t == NULL) return t;
    N.left = N.right = (Tree*)NULL;
    l = r = &N;
    N.weight = t->weight;
    
    y = t;
    for (;;) {
        if (i < y->item && (y->prev!=NULL && i<=y->prev->item)) {
            assert(i <= y->prev->item);
            if(y->right != NULL) right += y->right->weight;
            if(y->left == NULL) break;
            right += y->nodeWt;
            y = y->left;
        } else if (i > y->item) {
            if(y->left != NULL) left += y->left->weight;
            if(y->right == NULL) break;
            left += y->nodeWt;
            y = y->right;
        }
        else {
            /* i is within the block of y */
            if(y->right != NULL) right += y->right->weight;
            if(y->left != NULL) left += y->left->weight;
            break;
        }
    }
    for (;;) {
        if (i < t->item && (t->prev!=NULL && i<=t->prev->item)) {
            if (t->left == NULL) break;
            if (i < t->left->item && (t->left->prev!=NULL && i<=t->left->prev->item)) {
                y = t->left;     /* rotate right */
                t->left = y->right;
                y->right = t;
                /* t->weight--; */
                t->weight -= y->nodeWt;
                t = y;
                if (t->left == NULL) break;
                t->right->weight -= t->left->weight;
            }
            t->weight = right;
            /* right--; */
            right -= t->nodeWt;
            if(t->right != NULL) right -= t->right->weight;
            r->left = t;                               /* link right */
            r = t;
            t = t->left;
        } else if (i > t->item) {
            if (t->right == NULL) break;
            if (i > t->right->item) {
                y = t->right;                          /* rotate left */
                t->right = y->left;
                y->left = t;
                /* t->weight--; */
                t->weight -= y->nodeWt;
                t = y;
                if (t->right == NULL) break;
                t->left->weight -= t->right->weight;
            }
            t->weight = left;
            /* left--; */
            left -= t->nodeWt;
            if(t->left != NULL) left -= t->left->weight;
            l->right = t;                              /* link left */
            l = t;
            t = t->right;
        } else {
            break;
        }
    }
    l->right = t->left;                                /* assemble */
    r->left = t->right;
    t->left = N.right;
    t->right = N.left;
    t->weight = N.weight;
    
    *dis = t->nodeWt/2;
    if (t->right!=NULL) *dis += t->right->weight;
    
    return t;
}


Tree *ScaleTree::ScaleTreeInsertAtFront(unsigned long long blockEnd, Tree * t, Tree *newNode) {
    Tree * newOne, * prev  /*, * next*/;
    //unsigned useless;
    if (newNode==NULL) {
        //newOne = (Tree *) malloc (sizeof (Tree));
        newOne = MallocNode(); // by Zhangchl
        if (newOne == NULL) {
            printf("Ran out of space\n");
            exit(1);
        }
        sizeTrace ++;
    }
    else newOne = newNode;
    
    newOne->item = blockEnd;
    newOne->nodeWt = 1;
    newOne->maxSize = 1;
    if (t == NULL) {
        newOne->left = newOne->right = (Tree*)NULL;
        newOne->weight = 1;
        newOne->prev = newOne->next = (Tree*)NULL;
        return newOne;
    }
    
    /* Insert at the front of the tree */
    assert(blockEnd > t->item);
    newOne->weight = t->weight + 1;
    newOne->left = t;
    newOne->right = (Tree*)NULL;
    
    /* find prev and next */
    newOne->next = (Tree*)NULL;
    prev = newOne->left;
    if (prev!=NULL)
        while (prev->right!=NULL) prev = prev->right;
    newOne->prev = prev;
    if (prev!=NULL) prev->next = newOne;
    
    /* printf("insert: new %d, prev %d, next %d\n", new->item,
     new->prev!=NULL?new->prev->item: -1,
     new->next!=NULL? new->next->item: -1); */
    
    return newOne;
}

Tree *ScaleTree::QueryScaleTree(unsigned long long oldCyc, unsigned long long newCyc, Tree *t, unsigned *dis) {
    unsigned useless;
    Tree *tmp, *recycleNode = (Tree*)NULL;
    unsigned rightChildWt;
    if (oldCyc == newCyc) {
        t = ScaleTreeInsertAtFront(newCyc, t, (Tree*)NULL);
        return t;
    }
    t = ScaleTreeSplay(oldCyc, t, dis);
    assert(oldCyc <= t->item);
    if (t->prev!=NULL && oldCyc <= t->prev->item)
        assert(0);
    /* set the size of t */
    if (t->right!=NULL) rightChildWt = t->right->weight;
    else rightChildWt = 0;
    
    t->maxSize = (int) (rightChildWt * errorRate);
    if (t->maxSize==0) t->maxSize = 1;
    
    /* delete the oldCyc, merge nodes if necessary */
    t->nodeWt --;  t->weight --;
    assert(t->nodeWt>=0);
    if (t->nodeWt <= (t->maxSize >> 1)) {
        if (t->prev!=NULL && t->prev->nodeWt + t->nodeWt <= t->maxSize) {
            t->left = ScaleTreeSplay(t->prev->item, t->left, &useless);
            assert(t->left->right==NULL);  /* largest of the left subtree */
            assert(t->left==t->prev);
            t->left->right = t->right;  /* new tree */
            if (t->nodeWt > 0)   /* otherwise, t is empty */
                t->left->item = t->item;    /* merge  */
            t->left->nodeWt += t->nodeWt;
            t->left->weight += t->nodeWt;
            if (t->right!=NULL) t->left->weight += t->right->weight;
            t->left->next = t->next;    /* new neighbors */
            if (t->next!=NULL)
                t->next->prev = t->left;
            tmp = t;
            t = t->left;
            if (recycleNode == NULL) recycleNode = tmp;
            else {
                //free(tmp);
                FreeNode(tmp);
                sizeTrace --;
            }
        }
        if (t->prev!=NULL) {
            t->prev->maxSize = (int) ((rightChildWt + t->nodeWt) * errorRate);
            if (t->prev->maxSize==0) t->prev->maxSize = 1;
        }
        if (t->next!=NULL) {
            t->next->maxSize = (int) ((rightChildWt - t->next->nodeWt)
                                      * errorRate);
            if (t->next->maxSize==0) t->next->maxSize = 1;
        }
        if (t->next!=NULL && t->next->nodeWt+t->nodeWt <= t->next->maxSize) {
            /* merge next with me */
            t->right = ScaleTreeSplay(t->next->item, t->right, &useless);
            assert(t->right->left==NULL);
            assert(t->right==t->next);
            t->right->left = t->left;   /* new tree */
            t->right->nodeWt += t->nodeWt; /* merge    */
            t->right->weight += t->nodeWt;
            if (t->left!=NULL) t->right->weight += t->left->weight;
            t->right->prev = t->prev;    /* new neighbors */
            if (t->prev!=NULL)
                t->prev->next = t->right;
            tmp = t;
            t = t->right;
            if (recycleNode == NULL) recycleNode = tmp;
            else {
                //free(tmp);
                FreeNode(tmp);
                sizeTrace --;
            }
        }
    }
    
    /* ytzhong: only one address has been accessed and re-accessed,
     the nodeWt could be zero after compaction
     */
    /*  if (t->nodeWt == 0)
     assert(0);
     */
    
    
    /* insert newCyc */
    t = ScaleTreeInsertAtFront(newCyc, t, recycleNode);
    
    return t;
}

Tree *ScaleTree::CompactScaleTree(Tree *root) {
    if (root==NULL)
        return root;
    Tree *cur, *prev;
    unsigned priorWt = 0;
    unsigned totalWt = root->weight;
    /* from the end of trace forward, we set each node with correct
     maxSize and merge with its left neighbor if possible */
    assert(root->right==NULL);  /* root is the most recent */
    cur = root;
    while (cur!=NULL) {
        cur->right = (Tree*)NULL;  /* make it a list */
        cur->weight = totalWt - priorWt;
        cur->maxSize = (int) (priorWt * errorRate);
        if (cur->maxSize==0) cur->maxSize = 1;
        if (cur->prev!=NULL && cur->prev->nodeWt + cur->nodeWt <= cur->maxSize) {
            prev = cur->prev;
            cur->nodeWt += prev->nodeWt;
            cur->prev = prev->prev;
            if (cur->prev!=NULL)
                cur->prev->next = cur;
            cur->left = cur->prev;
            //free(prev);
            FreeNode(prev);
            sizeTrace --;
        }
        else {
            cur->left = cur->prev;
            priorWt += cur->nodeWt;
            cur = cur->left;
        }
    }
    return root;
}
