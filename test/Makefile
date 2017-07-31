#===- runtime/Makefile -------------------------*- Makefile -*-===##

# Relative path to the top of the source tree.
#
LEVEL=..

#LLVM_BIN_DIR=/u/cs255/build-llvm-38/bin
CLANG_BIN_DIR=/Users/dongchen/tools/llvm_xcode/Debug/bin
LLVM_BIN_DIR=/Users/dongchen/tools/llvm_xcode/Debug/bin

OPT= $(LLVM_BIN_DIR)/opt
LLC= $(LLVM_BIN_DIR)/llc
CC= $(CLANG_BIN_DIR)/clang

check: ./src/stencil.c
	$(CC) -O0 -c -emit-llvm ./src/stencil.c -o ./src_bc/stencil.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <./src_bc/stencil.bc> ./src_bc/stencil.bc.opt

correlation: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/correlation.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/correlation.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/correlation.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/correlation.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/correlation.bc.opt

convolution_3d: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_3d.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_3d.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_3d.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_3d.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_3d.bc.opt

convolution_2d: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_2d.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_2d.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_2d.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_2d.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/convolution_2d.bc.opt

syrk: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/syrk.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/syrk.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/syrk.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/syrk.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/syrk.bc.opt

mvt: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/mvt.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/mvt.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/mvt.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/mvt.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/mvt.bc.opt

gesummv: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gesummv.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gesummv.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gesummv.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gesummv.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gesummv.bc.opt

gemver: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemver.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemver.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemver.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemver.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemver.bc.opt

gemm: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemm.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemm.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemm.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemm.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/gemm.bc.opt

bicg: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/bicg.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/bicg.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/bicg.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/bicg.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/bicg.bc.opt

atax: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/atax.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/atax.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/atax.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/atax.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/atax.bc.opt

3mm: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/3mm.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/3mm.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/3mm.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/3mm.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/3mm.bc.opt

covariance: /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/covariance.c
	$(CC) -O0 -c -emit-llvm /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/covariance.c -o /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/covariance.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP </Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/covariance.bc> /Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/covariance.bc.opt

clean:
	rm -f *.bc *.opt *.o test

