
 /* Start to analysis array index
Array index info
mean.addr j
data.addr ((i * 1024) + j)
mean.addr j
mean.addr j
mean.addr j
mean.addr j
stddev.addr j
data.addr ((i * 1024) + j)
mean.addr j
data.addr ((i * 1024) + j)
mean.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
mean.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
stddev.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
symmat.addr ((j1 * 1024) + j1)
symmat.addr ((j1 * 1024) + j2)
data.addr ((i * 1024) + j1)
data.addr ((i * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j2 * 1024) + j1)
symmat.addr 1048575

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %m
i32 %n
double* %data
double* %mean
double* %stddev
double* %symmat

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
----array access   %14 = load double, double* %arrayidx10, align 8
----array access   store double %div, double* %arrayidx10, align 8
--j
--Loop Bound: (0, 1024)
----array access   store double 0.000000e+00, double* %arrayidx18, align 8
----i
----Loop Bound: (0, 1024)
------array access   %23 = load double, double* %arrayidx25, align 8
------array access   %26 = load double, double* %arrayidx27, align 8
------array access   %30 = load double, double* %arrayidx31, align 8
------array access   %33 = load double, double* %arrayidx33, align 8
------array access   %36 = load double, double* %arrayidx37, align 8
------array access   store double %add38, double* %arrayidx37, align 8
----array access   %40 = load double, double* %arrayidx43, align 8
----array access   store double %div44, double* %arrayidx43, align 8
----array access   store double %conv, double* %arrayidx46, align 8
----array access   %47 = load double, double* %arrayidx48, align 8
----array access   %50 = load double, double* %arrayidx52, align 8
----array access   store double %cond, double* %arrayidx54, align 8
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access   %58 = load double, double* %arrayidx67, align 8
------array access   %62 = load double, double* %arrayidx71, align 8
------array access   store double %sub72, double* %arrayidx71, align 8
------array access   %65 = load double, double* %arrayidx75, align 8
------array access   %69 = load double, double* %arrayidx80, align 8
------array access   store double %div81, double* %arrayidx80, align 8
--j1
--Loop Bound: (0, 1024)
----array access   store double 1.000000e+00, double* %arrayidx95, align 8
----j2
----Loop Bound: ((j1 + 1), 1024)
------array access   store double 0.000000e+00, double* %arrayidx104, align 8
------i
------Loop Bound: (0, 1024)
--------array access   %85 = load double, double* %arrayidx112, align 8
--------array access   %89 = load double, double* %arrayidx116, align 8
--------array access   %93 = load double, double* %arrayidx121, align 8
--------array access   store double %add122, double* %arrayidx121, align 8
------array access   %98 = load double, double* %arrayidx129, align 8
------array access   store double %98, double* %arrayidx133, align 8
--array access   store double 1.000000e+00, double* %arrayidx140, align 8

Finish analysis loops */ 
 // Start to generating Static Sampling Code
Loop at depth 2 containing: %for.cond1<header><exiting>,%for.body3,%for.inc<latch>
Loop at depth 2 containing: %for.cond19<header><exiting>,%for.body21,%for.inc39<latch>
Loop at depth 2 containing: %for.cond62<header><exiting>,%for.body65,%for.inc82<latch>
Loop at depth 2 containing: %for.cond97<header><exiting>,%for.body100,%for.cond105,%for.end125,%for.inc134<latch>,%for.body108,%for.inc123
    Loop at depth 3 containing: %for.cond105<header><exiting>,%for.body108,%for.inc123<latch>
libc++abi.dylib: terminating with uncaught exception of type std::invalid_argument: stoi: no conversion
0  opt                      0x000000010b288bbc llvm::sys::PrintStackTrace(llvm::raw_ostream&) + 60
1  opt                      0x000000010b289139 PrintStackTraceSignalHandler(void*) + 25
2  opt                      0x000000010b285269 llvm::sys::RunSignalHandlers() + 425
3  opt                      0x000000010b289802 SignalHandler(int) + 354
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
14 LLVMSymFP.dylib          0x000000011085092f ssCodeGen::StaticSamplingCodeGen::initRefCntOfLoop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*) + 735
15 LLVMSymFP.dylib          0x0000000110850c64 ssCodeGen::StaticSamplingCodeGen::initRefCntOfLoop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*) + 1556
16 LLVMSymFP.dylib          0x0000000110853469 ssCodeGen::StaticSamplingCodeGen::init(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*) + 73
17 LLVMSymFP.dylib          0x000000011085358c ssCodeGen::StaticSamplingCodeGen::runOnFunction(llvm::Function&) + 236
18 opt                      0x000000010a93badf llvm::FPPassManager::runOnFunction(llvm::Function&) + 399
19 opt                      0x000000010a93bfe5 llvm::FPPassManager::runOnModule(llvm::Module&) + 117
20 opt                      0x000000010a93cdb4 (anonymous namespace)::MPPassManager::runOnModule(llvm::Module&) + 2196
21 opt                      0x000000010a93c2a6 llvm::legacy::PassManagerImpl::run(llvm::Module&) + 342
22 opt                      0x000000010a93db01 llvm::legacy::PassManager::run(llvm::Module&) + 33
23 opt                      0x000000010895d475 main + 25429
24 libdyld.dylib            0x00007fffd9108235 start + 1
25 libdyld.dylib            0x0000000000000004 start + 653229520
Stack dump:
0.	Program arguments: /Users/dongchen/tools/llvm_xcode/Debug/bin/opt -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP 
1.	Running pass 'Function Pass Manager' on module '<stdin>'.
2.	Running pass 'static sampling code generating pass' on function '@correlation'
