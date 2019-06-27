/*BEGIN_LEGAL 
Intel Open Source License 

Copyright (c) 2002-2009 Intel Corporation. All rights reserved.
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.  Redistributions
in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.  Neither the name of
the Intel Corporation nor the names of its contributors may be used to
endorse or promote products derived from this software without
specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE INTEL OR
ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
END_LEGAL */
/* ===================================================================== */
/*
  @ORIGINAL_AUTHOR: Robert Muth
*/

/* ===================================================================== */
/*! @file
 *  This file contains an ISA-portable PIN tool for counting dynamic instructions
 */

//#include "pin.H"
#include <iostream>
//#include <fstream>
#include "scaletree-spatial.h"
#include "stdlib.h"

/* ===================================================================== */
/* For Reuse Distance Analysis */
/* ===================================================================== */
//std::ofstream ResultFile;
unsigned long long matrix[NumCounters][NumCounters];
ScaleTree *tree1 = (ScaleTree*)NULL;
ScaleTree *tree2 = (ScaleTree*)NULL;

void ScaleTreeInit(ScaleTree *tree, int allocBuf) {
    tree->trace = (Tree*)NULL;
    tree->numInstr = 0;
    tree->logblocksize = 0;
    tree->buffersize = BUFFERSIZE;
    tree->logTrace = 0;
    tree->instr = 0;
    tree->errorRate = 0.001;
    if (allocBuf != 0) {
	tree->buffer = (unsigned long*)malloc(tree->buffersize*sizeof(unsigned long));
	tree->rd2 = (unsigned*)malloc(tree->buffersize*sizeof(unsigned));
    } else {
	tree->buffer = (unsigned long*)NULL;
	tree->rd2 = (unsigned*)NULL;
    }

    tree->buffer = (unsigned long*)malloc(tree->buffersize*sizeof(unsigned long));
    tree->hashUsed = 0;
    tree->hashallocated = 0;
    tree->numData = 0;
    tree->curCycle = 0;
    tree->power = 1;
    tree->bufferindex = 0;
    tree->buffersprocessed = 0;
    tree->sizeTrace = 0;
    tree->freeNodeList = (Tree*)NULL;

    tree->Init();
}

int get_bin_number(unsigned rd) {
    int val = 0;
    if (rd == (unsigned)-1) {
	return NumCounters-1;
    }
    while ((rd > 0) && (val < NumCounters-2)) {
	val++;
	rd=rd>>1;
    }
    return val;
}

void processBuffer() {
    //ResultFile << "processBuffer() is being called\t" << tree1->bufferindex << endl;
    tree1->buffersprocessed++;

    if (tree1->buffersprocessed%128 == 0) {
	printf("processing %ld buffers, bufferindex=%d\n",tree1->buffersprocessed, tree1->bufferindex);
	tree1->PrintSize();
    }

    if(tree1->logTrace) {
	for(int i=0;i<tree1->bufferindex;i++) {
	    printf("%p\n",(void*)(tree1->buffer[i] >> tree1->logblocksize));
	}
    }

    unsigned rd1,rd2;

    for (int i=0;i<tree2->bufferindex;i++) {
	rd2 = tree2->DataAccess(tree2->buffer[i]>>tree2->logblocksize);
	tree1->rd2[i] = rd2;
    }

    for (int i=0;i<tree1->bufferindex;i++) {
	int bin1,bin2;
	rd1 = tree1->DataAccess(tree1->buffer[i]>>tree1->logblocksize);
	bin1 = get_bin_number(rd1);
	bin2 = get_bin_number(tree1->rd2[i]);
	matrix[bin1][bin2]++;
    }
    tree1->bufferindex=0;
    tree2->bufferindex=0;
}


/* ===================================================================== */
/* Commandline Switches */
/* ===================================================================== */

//KNOB<string> KnobResultFile(KNOB_MODE_WRITEONCE, "pintool",
//			    "o", "reda-spatial-results.txt", "specify result file name");
//KNOB<string> KnobTraceFile(KNOB_MODE_WRITEONCE, "pintool",
//			    "o", "reda-spatial-trace.txt", "specify trace file name");

/* ===================================================================== */

//INT32 Usage() {
//    cerr << "This tool prints out the number of dynamic instructions executed to stderr.\n";

//    cerr << KNOB_BASE::StringKnobSummary();

//    cerr << endl;

//    return -1;
//}

/* ===================================================================== */

//VOID docount() {
//    ins_count++;
//}

// Print a memory read record
//VOID RecordMemRead(VOID * ip, VOID * addr) {
//    tree2->bufferindex++;
//    tree1->buffer[tree1->bufferindex++] = (unsigned long)addr;
//    if (tree1->bufferindex == tree1->buffersize)
//	processBuffer();
//    //    read_count++;
//    //TraceFile << ip << ": R " << addr << endl;
//}

// Print a memory write record
//VOID RecordMemWrite(VOID * ip, VOID * addr) {
//    tree2->bufferindex++;
//    tree1->buffer[tree1->bufferindex++] = (unsigned long)addr;
//    if (tree1->bufferindex == tree1->buffersize)
//	processBuffer();
//   //    write_count++;
//    //TraceFile << ip << ": W " << addr << endl;
//}

void rtTmpAccess(int addr) {

	addr = addr * 8;

    tree2->bufferindex++;
	tree1->buffer[tree1->bufferindex++] = (unsigned long)addr;
	if (tree1->bufferindex == tree1->buffersize)
		processBuffer();
}


/* ===================================================================== */

//VOID Instruction(INS ins, VOID *v) {
    //INS_InsertCall(ins, IPOINT_BEFORE, (AFUNPTR)docount, IARG_END);

    // Instruments loads using a predicated call, i.e.
    // the call happens iff the load will be actually executed.
    // The IA-64 architecture has explicitly predicated instructions. 
    // On the IA-32 and Intel(R) 64 architectures conditional moves and REP 
    // prefixed instructions appear as predicated instructions in Pin.
//    if (INS_IsMemoryRead(ins)) {
//	INS_InsertPredicatedCall(ins, IPOINT_BEFORE, (AFUNPTR)RecordMemRead,
//				 IARG_INST_PTR, IARG_MEMORYREAD_EA, 
//				 IARG_END);
//    }

    // Instruments stores using a predicated call, i.e.
    // the call happens iff the store will be actually executed
//    if (INS_IsMemoryWrite(ins))	{
//	INS_InsertPredicatedCall(ins, IPOINT_BEFORE, (AFUNPTR)RecordMemWrite,
//				 IARG_INST_PTR, IARG_MEMORYWRITE_EA, IARG_END);
//    }
//}

/* ===================================================================== */

//VOID Fini(INT32 code, VOID *v) {
void FiniRD() {
  //cerr <<  "Count " << ins_count  << endl;
  //ResultFile << "Inst count: " << ins_count << endl;
  //ResultFile << "MemRead count: " << read_count << endl;
  //ResultFile << "MemWrite count: " << write_count << endl;
    processBuffer();

    printf("Starting to print results.\n");
    for (int i=0;i<=NumCounters;i++) {
	unsigned long long sum = 0;
	unsigned long long effect_sum = 0;
	for (int j=0;j<=i;j++) {
	    sum += matrix[i][j];
	    if (j <= i-REDUCED_THRESHOLD)
		effect_sum += matrix[i][j];
	}
	if (sum != 0)
	    tree1->slq[i] = effect_sum*2.0/sum;
    }
//    tree1->PrintResults(&ResultFile);
	tree1->PrintResults();
//	tree2->PrintResults();

//    ResultFile.close();
  //TraceFile.close();
}

/* ===================================================================== */

//int main(int argc, char *argv[]) {
int InitRD() {
//    if( PIN_Init(argc,argv) ) {
//        return Usage();
//    }

    //initialize the matrix for the effect of doubling block size
    int i,j;
    for (i=0;i<NumCounters;i++)
	for (j=0;j<NumCounters;j++)
	    matrix[i][j]=0;

    if(tree1 == NULL) {
		tree1 = (ScaleTree*)malloc(sizeof(ScaleTree));
		ScaleTreeInit(tree1,1);
		tree1->logblocksize = LOG_BLOCK_SIZE;
    }

    if (tree2 == NULL) {
		tree2 = (ScaleTree*)malloc(sizeof(ScaleTree));
		ScaleTreeInit(tree2,0);
		tree1->logblocksize = LOG_BLOCK_SIZE+1;
		tree2->buffer = tree1->buffer;
    }
    
//    ResultFile.open(KnobResultFile.Value().c_str());
    //TraceFile.open(KnobTraceFile.Value().c_str());

//    INS_AddInstrumentFunction(Instruction, 0);
//    PIN_AddFiniFunction(Fini, 0);

    // Never returns
//    PIN_StartProgram();
  
    return 0;
}


/* ===================================================================== */
/* eof */
/* ===================================================================== */
