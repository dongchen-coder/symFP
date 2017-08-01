#===- runtime/Makefile -------------------------*- Makefile -*-===##

# Relative path to the top of the source tree.
#
LEVEL=..

#LLVM_BIN_DIR=/u/cs255/build-llvm-38/bin
CLANG_BIN_DIR=/Users/dongchen/tools/llvm_xcode/Debug/bin
LLVM_BIN_DIR=/Users/dongchen/tools/llvm_xcode/Debug/bin

TEST_DIR=/Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/polyBench


OPT= $(LLVM_BIN_DIR)/opt
LLC= $(LLVM_BIN_DIR)/llc
CC= $(CLANG_BIN_DIR)/clang

polyBench: check correlation convolution_3d convolution_2d syrk mvt gesummv gemver gemm bicg atax covariance

check: $(TEST_DIR)/stencil.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/stencil.c -o $(TEST_DIR)/stencil.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/stencil.bc> $(TEST_DIR)/stencil.bc.opt

correlation: $(TEST_DIR)/correlation.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/correlation.c -o $(TEST_DIR)/correlation.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/correlation.bc> $(TEST_DIR)/correlation.bc.opt

convolution_3d: $(TEST_DIR)/convolution_3d.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/convolution_3d.c -o $(TEST_DIR)/convolution_3d.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/convolution_3d.bc> $(TEST_DIR)/convolution_3d.bc.opt

convolution_2d: $(TEST_DIR)/convolution_2d.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/convolution_2d.c -o $(TEST_DIR)/convolution_2d.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/convolution_2d.bc> $(TEST_DIR)/convolution_2d.bc.opt

syrk: $(TEST_DIR)/syrk.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/syrk.c -o $(TEST_DIR)/syrk.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/syrk.bc> $(TEST_DIR)/syrk.bc.opt

mvt: $(TEST_DIR)/mvt.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/mvt.c -o $(TEST_DIR)/mvt.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/mvt.bc> $(TEST_DIR)/mvt.bc.opt

gesummv: $(TEST_DIR)/gesummv.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/gesummv.c -o $(TEST_DIR)/gesummv.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/gesummv.bc> $(TEST_DIR)/gesummv.bc.opt

gemver: $(TEST_DIR)/gemver.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/gemver.c -o $(TEST_DIR)/gemver.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/gemver.bc> $(TEST_DIR)/gemver.bc.opt

gemm: $(TEST_DIR)/gemm.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/gemm.c -o $(TEST_DIR)/gemm.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/gemm.bc> $(TEST_DIR)/gemm.bc.opt

bicg: $(TEST_DIR)/bicg.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/bicg.c -o $(TEST_DIR)/bicg.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/bicg.bc> $(TEST_DIR)/bicg.bc.opt

atax: $(TEST_DIR)/atax.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/atax.c -o $(TEST_DIR)/atax.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/atax.bc> $(TEST_DIR)/atax.bc.opt

3mm: $(TEST_DIR)/3mm.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/3mm.c -o $(TEST_DIR)/3mm.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/3mm.bc> $(TEST_DIR)/3mm.bc.opt

covariance: $(TEST_DIR)/covariance.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/covariance.c -o $(TEST_DIR)/covariance.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/covariance.bc> $(TEST_DIR)/covariance.bc.opt

clean:
	rm -f *.bc *.opt *.o test

