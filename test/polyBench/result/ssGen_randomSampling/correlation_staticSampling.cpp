
 /* Start to analysis array index
0  opt                      0x0000000107d12bbc llvm::sys::PrintStackTrace(llvm::raw_ostream&) + 60
1  opt                      0x0000000107d13139 PrintStackTraceSignalHandler(void*) + 25
2  opt                      0x0000000107d0f269 llvm::sys::RunSignalHandlers() + 425
3  opt                      0x0000000107d13802 SignalHandler(int) + 354
4  libsystem_platform.dylib 0x00007fffd9317b3a _sigtramp + 26
5  libsystem_platform.dylib 0x00007fff5a840468 _sigtramp + 2169669960
6  LLVMSymFP.dylib          0x000000010d263af0 llvm::UnaryInstruction::getOperand(unsigned int) const + 160
7  LLVMSymFP.dylib          0x000000010d25b34a idxAnalysis::IndexAnalysis::computeExpression(llvm::Instruction*) + 2890
8  LLVMSymFP.dylib          0x000000010d265b08 idxAnalysis::IndexAnalysis::findAllArrayAccesses(llvm::Function&) + 4040
9  LLVMSymFP.dylib          0x000000010d26725d idxAnalysis::IndexAnalysis::runOnFunction(llvm::Function&) + 61
10 opt                      0x00000001073c5adf llvm::FPPassManager::runOnFunction(llvm::Function&) + 399
11 opt                      0x00000001073c5fe5 llvm::FPPassManager::runOnModule(llvm::Module&) + 117
12 opt                      0x00000001073c6db4 (anonymous namespace)::MPPassManager::runOnModule(llvm::Module&) + 2196
13 opt                      0x00000001073c62a6 llvm::legacy::PassManagerImpl::run(llvm::Module&) + 342
14 opt                      0x00000001073c7b01 llvm::legacy::PassManager::run(llvm::Module&) + 33
15 opt                      0x00000001053e7475 main + 25429
16 libdyld.dylib            0x00007fffd9108235 start + 1
17 libdyld.dylib            0x0000000000000004 start + 653229520
Stack dump:
0.	Program arguments: /Users/dongchen/tools/llvm_xcode/Debug/bin/opt -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP 
1.	Running pass 'Function Pass Manager' on module '<stdin>'.
2.	Running pass 'Array index analysis Pass' on function '@correlation'
