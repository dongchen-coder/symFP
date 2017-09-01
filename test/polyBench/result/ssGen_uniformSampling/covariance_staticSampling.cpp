
 /* Start to analysis array index
Array index info
data.addr ((i * 1024) + j)
mean.addr j
mean.addr j
mean.addr j
mean.addr j
mean.addr j
mean.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
symmat.addr ((j1 * 1024) + j2)
data.addr ((i * 1024) + j1)
data.addr ((i * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j2 * 1024) + j1)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %m
i32 %n
double* %data
double* %symmat
double* %mean

 Finish to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--j
--Loop Bound: (0, 1024)
----array access   store double 0.000000e+00, double* %arrayidx, align 8
----i
----Loop Bound: (0, 1024)
------array access   %7 = load double, double* %arrayidx5, align 8
------array access   %10 = load double, double* %arrayidx7, align 8
------array access   store double %add8, double* %arrayidx7, align 8
----array access   %15 = load double, double* %arrayidx10, align 8
----array access   store double %div, double* %arrayidx10, align 8
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access   %21 = load double, double* %arrayidx21, align 8
------array access   %25 = load double, double* %arrayidx25, align 8
------array access   store double %sub, double* %arrayidx25, align 8
--j1
--Loop Bound: (0, 1024)
----j2
----Loop Bound: (j1, 1024)
------array access   store double 0.000000e+00, double* %arrayidx41, align 8
------i
------Loop Bound: (0, 1024)
--------array access   %38 = load double, double* %arrayidx48, align 8
--------array access   %42 = load double, double* %arrayidx52, align 8
--------array access   %46 = load double, double* %arrayidx57, align 8
--------array access   store double %add58, double* %arrayidx57, align 8
------array access   %51 = load double, double* %arrayidx65, align 8
------array access   store double %51, double* %arrayidx69, align 8

Finish analysis loops */ 
 // Start to generating Static Sampling Code
libc++abi.dylib: terminating with uncaught exception of type std::invalid_argument: stoi: no conversion
0  opt                      0x0000000103e53bbc llvm::sys::PrintStackTrace(llvm::raw_ostream&) + 60
1  opt                      0x0000000103e54139 PrintStackTraceSignalHandler(void*) + 25
2  opt                      0x0000000103e50269 llvm::sys::RunSignalHandlers() + 425
3  opt                      0x0000000103e54802 SignalHandler(int) + 354
4  libsystem_platform.dylib 0x00007fffd9317b3a _sigtramp + 26
5  libsystem_platform.dylib 0x00000004e2000240 _sigtramp + 147752736
6  libsystem_c.dylib        0x00007fffd919c420 abort + 129
7  libc++abi.dylib          0x00007fffd7cef94a __cxa_bad_cast + 0
8  libc++abi.dylib          0x00007fffd7d14c17 default_terminate_handler() + 243
9  libobjc.A.dylib          0x00007fffd8824713 _objc_terminate() + 124
10 libc++abi.dylib          0x00007fffd7d11d49 std::__terminate(void (*)()) + 8
11 libc++abi.dylib          0x00007fffd7d117be __cxxabiv1::exception_cleanup_func(_Unwind_Reason_Code, _Unwind_Exception*) + 0
12 libc++.1.dylib           0x00007fffd7cde6d0 long std::__1::(anonymous namespace)::as_integer_helper<long, std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> >, long (*)(char const*, char**, int)>(std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > const&, std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > const&, unsigned long*, int, long (*)(char const*, char**, int)) + 297
13 libc++.1.dylib           0x00007fffd7cdbc40 std::__1::stoi(std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > const&, unsigned long*, int) + 80
14 LLVMSymFP.dylib          0x000000010941b941 ssCodeGen::StaticSamplingCodeGen::initRefCntOfLoop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*) + 705
15 LLVMSymFP.dylib          0x000000010941bc70 ssCodeGen::StaticSamplingCodeGen::initRefCntOfLoop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*) + 1520
16 LLVMSymFP.dylib          0x000000010941e479 ssCodeGen::StaticSamplingCodeGen::init(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*) + 73
17 LLVMSymFP.dylib          0x000000010941e59c ssCodeGen::StaticSamplingCodeGen::runOnFunction(llvm::Function&) + 236
18 opt                      0x0000000103506adf llvm::FPPassManager::runOnFunction(llvm::Function&) + 399
19 opt                      0x0000000103506fe5 llvm::FPPassManager::runOnModule(llvm::Module&) + 117
20 opt                      0x0000000103507db4 (anonymous namespace)::MPPassManager::runOnModule(llvm::Module&) + 2196
21 opt                      0x00000001035072a6 llvm::legacy::PassManagerImpl::run(llvm::Module&) + 342
22 opt                      0x0000000103508b01 llvm::legacy::PassManager::run(llvm::Module&) + 33
23 opt                      0x0000000101528475 main + 25429
24 libdyld.dylib            0x00007fffd9108235 start + 1
25 libdyld.dylib            0x0000000000000004 start + 653229520
Stack dump:
0.	Program arguments: /Users/dongchen/tools/llvm_xcode/Debug/bin/opt -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP 
1.	Running pass 'Function Pass Manager' on module '<stdin>'.
2.	Running pass 'static sampling code generating pass' on function '@covariance'
