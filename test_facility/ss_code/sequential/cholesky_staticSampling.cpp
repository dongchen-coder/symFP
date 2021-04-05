
 /* Start to analysis array index
Array index info: Total number of references: 13
A.addr ((i * 1024) + k)
A.addr ((j * 1024) + k)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + i)
A.addr ((j * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + i)
A.addr ((i * 1024) + i)
A.addr ((i * 1024) + i)
BC Array cost info: Total number of arrays: 1
A.addr 24
BC Average cost: 2.400000e+01

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, i)
----Loop inc: (j + 1)
----Loop predicate: <
------k
------Loop Bound: (0, j)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 1024) + k)
--------array access A.addr ((j * 1024) + k)
--------array access A.addr ((i * 1024) + j)
--------array access A.addr ((i * 1024) + j)
------array access A.addr ((j * 1024) + j)
------array access A.addr ((i * 1024) + j)
------array access A.addr ((i * 1024) + j)
----k
----Loop Bound: (0, i)
----Loop inc: (k + 1)
----Loop predicate: <
------array access A.addr ((i * 1024) + k)
------array access A.addr ((i * 1024) + k)
------array access A.addr ((i * 1024) + i)
------array access A.addr ((i * 1024) + i)
----array access A.addr ((i * 1024) + i)
----array access A.addr ((i * 1024) + i)

Finish analysis loops */ 
terminate called after throwing an instance of 'std::invalid_argument'
  what():  stoi
PLEASE submit a bug report to https://bugs.llvm.org/ and include the crash backtrace.
Stack dump:
0.	Program arguments: /localdisk/fliu14/llvm/fzllvm/build/bin/opt -load ../build/lib/libLLVMSPS.so -sps -spsrate=0.02 
1.	Running pass 'Function Pass Manager' on module '<stdin>'.
2.	Running pass 'loop tree transform for parallel program Pass' on function '@cholesky'
 #0 0x0000000001b0d1fa llvm::sys::PrintStackTrace(llvm::raw_ostream&) (/localdisk/fliu14/llvm/fzllvm/build/bin/opt+0x1b0d1fa)
 #1 0x0000000001b0b2f4 llvm::sys::RunSignalHandlers() (/localdisk/fliu14/llvm/fzllvm/build/bin/opt+0x1b0b2f4)
 #2 0x0000000001b0b465 SignalHandler(int) (/localdisk/fliu14/llvm/fzllvm/build/bin/opt+0x1b0b465)
 #3 0x00007fb5564b8a90 __restore_rt (/lib64/libpthread.so.0+0x14a90)
 #4 0x00007fb555f697d5 raise (/lib64/libc.so.6+0x3c7d5)
 #5 0x00007fb555f52895 abort (/lib64/libc.so.6+0x25895)
 #6 0x00007fb5562f7941 (/lib64/libstdc++.so.6+0x9e941)
 #7 0x00007fb55630342c (/lib64/libstdc++.so.6+0xaa42c)
 #8 0x00007fb556303497 (/lib64/libstdc++.so.6+0xaa497)
 #9 0x00007fb556303749 (/lib64/libstdc++.so.6+0xaa749)
#10 0x00007fb5562fa350 std::__throw_invalid_argument(char const*) (/lib64/libstdc++.so.6+0xa1350)
#11 0x00007fb555dd3845 int __gnu_cxx::__stoa<long, int, char, int>(long (*)(char const*, char**, int), char const*, char const*, unsigned long*, int) (../build/lib/libLLVMSPS.so+0x150845)
#12 0x00007fb555dd3687 std::__cxx11::stoi(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, unsigned long*, int) (../build/lib/libLLVMSPS.so+0x150687)
#13 0x00007fb555dd2a83 loopTreeTransform::ParallelLoopTreeTransform::computeIterSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*) (.localalias) (../build/lib/libLLVMSPS.so+0x14fa83)
#14 0x00007fb555dd29a6 loopTreeTransform::ParallelLoopTreeTransform::computeIterSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*) (.localalias) (../build/lib/libLLVMSPS.so+0x14f9a6)
#15 0x00007fb555dd29a6 loopTreeTransform::ParallelLoopTreeTransform::computeIterSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*) (.localalias) (../build/lib/libLLVMSPS.so+0x14f9a6)
#16 0x00007fb555dd2d90 loopTreeTransform::ParallelLoopTreeTransform::computePerIterationSpace() (../build/lib/libLLVMSPS.so+0x14fd90)
#17 0x00007fb555dd352e loopTreeTransform::ParallelLoopTreeTransform::runOnFunction(llvm::Function&) (../build/lib/libLLVMSPS.so+0x15052e)
#18 0x00000000013010c3 llvm::FPPassManager::runOnFunction(llvm::Function&) (/localdisk/fliu14/llvm/fzllvm/build/bin/opt+0x13010c3)
#19 0x0000000001301711 llvm::FPPassManager::runOnModule(llvm::Module&) (/localdisk/fliu14/llvm/fzllvm/build/bin/opt+0x1301711)
#20 0x00000000013000ce llvm::legacy::PassManagerImpl::run(llvm::Module&) (/localdisk/fliu14/llvm/fzllvm/build/bin/opt+0x13000ce)
#21 0x00000000006a52c0 main (/localdisk/fliu14/llvm/fzllvm/build/bin/opt+0x6a52c0)
#22 0x00007fb555f54082 __libc_start_main (/lib64/libc.so.6+0x27082)
#23 0x000000000073240e _start (/localdisk/fliu14/llvm/fzllvm/build/bin/opt+0x73240e)
