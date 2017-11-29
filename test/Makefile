#===- runtime/Makefile -------------------------*- Makefile -*-===##

# Relative path to the top of the source tree.
#
LEVEL=..

#LLVM_BIN_DIR=/u/cs255/build-llvm-38/bin
CLANG_BIN_DIR=/Users/dongchen/tools/llvm_xcode/Debug/bin
LLVM_BIN_DIR=/Users/dongchen/tools/llvm_xcode/Debug/bin

TEST_DIR=/Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/test/polyBench
RESULT_DIR=$(TEST_DIR)/result_ref

OPT= $(LLVM_BIN_DIR)/opt
LLC= $(LLVM_BIN_DIR)/llc
CC= $(CLANG_BIN_DIR)/clang

polyBench: stencil convolution_2d syrk mvt gesummv gemver gemm bicg atax 3mm convolution_3d
# correlation covariance

time_codeGen: stencil_time_codeGen convolution_2d_time_codeGen syrk_time_codeGen mvt_time_codeGen gesummv_time_codeGen gemver_time_codeGen gemm_time_codeGen bicg_time_codeGen atax_time_codeGen 3mm_time_codeGen convolution_3d_time_codeGen

check: $(TEST_DIR)/stencil.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/stencil.c -o $(TEST_DIR)/stencil.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/stencil.bc> $(TEST_DIR)/stencil.bc.opt

stencil: $(TEST_DIR)/stencil.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/stencil.c -o $(TEST_DIR)/stencil.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/stencil.bc> $(TEST_DIR)/stencil.bc.opt 2> $(RESULT_DIR)/stencil_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/stencil_staticSampling.cpp -O3 -o $(RESULT_DIR)/stencil_staticSampling

stencil_time_codeGen: $(TEST_DIR)/stencil.c
	#echo "stencil" > $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/stencil.bc> $(TEST_DIR)/stencil.bc.opt 2> $(RESULT_DIR)/stencil_staticSampling.cpp ; }
# 2>> $(RESULT_DIR)/time.log 

correlation: $(TEST_DIR)/correlation.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/correlation.c -o $(TEST_DIR)/correlation.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/correlation.bc> $(TEST_DIR)/correlation.bc.opt 2> $(RESULT_DIR)/correlation_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/correlation_staticSampling.cpp -O3 -o $(RESULT_DIR)/correlation_staticSampling

correlation_time_codeGen: $(TEST_DIR)/correlation.c
	echo "correlation" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/correlation.bc> $(TEST_DIR)/correlation.bc.opt 2> $(RESULT_DIR)/correlation_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log 

convolution_3d: $(TEST_DIR)/convolution_3d.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/convolution_3d.c -o $(TEST_DIR)/convolution_3d.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/convolution_3d.bc> $(TEST_DIR)/convolution_3d.bc.opt 2> $(RESULT_DIR)/convolution_3d_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/convolution_3d_staticSampling.cpp -O3 -o $(RESULT_DIR)/convolution_3d_staticSampling

convolution_3d_time_codeGen: $(TEST_DIR)/convolution_3d.c
	echo "convolution_3d" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/convolution_3d.bc> $(TEST_DIR)/convolution_3d.bc.opt 2> $(RESULT_DIR)/convolution_3d_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log 

convolution_2d: $(TEST_DIR)/convolution_2d.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/convolution_2d.c -o $(TEST_DIR)/convolution_2d.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/convolution_2d.bc> $(TEST_DIR)/convolution_2d.bc.opt 2> $(RESULT_DIR)/convolution_2d_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/convolution_2d_staticSampling.cpp -O3 -o $(RESULT_DIR)/convolution_2d_staticSampling

convolution_2d_time_codeGen: $(TEST_DIR)/convolution_2d.c
	echo "convolution_2d" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/convolution_2d.bc> $(TEST_DIR)/convolution_2d.bc.opt 2> $(RESULT_DIR)/convolution_2d_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log 

syrk: $(TEST_DIR)/syrk.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/syrk.c -o $(TEST_DIR)/syrk.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/syrk.bc> $(TEST_DIR)/syrk.bc.opt 2> $(RESULT_DIR)/syrk_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/syrk_staticSampling.cpp -O3 -o $(RESULT_DIR)/syrk_staticSampling

syrk_time_codeGen: $(TEST_DIR)/syrk.c
	echo "syrk" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/syrk.bc> $(TEST_DIR)/syrk.bc.opt 2> $(RESULT_DIR)/syrk_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log

mvt: $(TEST_DIR)/mvt.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/mvt.c -o $(TEST_DIR)/mvt.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/mvt.bc> $(TEST_DIR)/mvt.bc.opt 2> $(RESULT_DIR)/mvt_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/mvt_staticSampling.cpp -O3 -o $(RESULT_DIR)/mvt_staticSampling

mvt_time_codeGen: $(TEST_DIR)/mvt.c
	echo "mvt" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/mvt.bc> $(TEST_DIR)/mvt.bc.opt 2> $(RESULT_DIR)/mvt_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log

gesummv: $(TEST_DIR)/gesummv.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/gesummv.c -o $(TEST_DIR)/gesummv.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/gesummv.bc> $(TEST_DIR)/gesummv.bc.opt 2> $(RESULT_DIR)/gesummv_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/gesummv_staticSampling.cpp -O3 -o $(RESULT_DIR)/gesummv_staticSampling
	
gesummv_time_codeGen: $(TEST_DIR)/gesummv.c
	echo "gesumv" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/gesummv.bc> $(TEST_DIR)/gesummv.bc.opt 2> $(RESULT_DIR)/gesummv_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log

gemver: $(TEST_DIR)/gemver.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/gemver.c -o $(TEST_DIR)/gemver.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/gemver.bc> $(TEST_DIR)/gemver.bc.opt 2> $(RESULT_DIR)/gemver_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/gemver_staticSampling.cpp -O3 -o $(RESULT_DIR)/gemver_staticSampling

gemver_time_codeGen: $(TEST_DIR)/gemver.c
	echo "gemver" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/gemver.bc> $(TEST_DIR)/gemver.bc.opt 2> $(RESULT_DIR)/gemver_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log

gemm: $(TEST_DIR)/gemm.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/gemm.c -o $(TEST_DIR)/gemm.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/gemm.bc> $(TEST_DIR)/gemm.bc.opt 2> $(RESULT_DIR)/gemm_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/gemm_staticSampling.cpp -O3 -o $(RESULT_DIR)/gemm_staticSampling

gemm_time_codeGen: $(TEST_DIR)/gemm.c
	echo "gemm" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/gemm.bc> $(TEST_DIR)/gemm.bc.opt 2> $(RESULT_DIR)/gemm_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log

bicg: $(TEST_DIR)/bicg.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/bicg.c -o $(TEST_DIR)/bicg.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/bicg.bc> $(TEST_DIR)/bicg.bc.opt 2> $(RESULT_DIR)/bicg_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/bicg_staticSampling.cpp -O3 -o $(RESULT_DIR)/bicg_staticSampling

bicg_time_codeGen: $(TEST_DIR)/bicg.c
	echo "bicg" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/bicg.bc> $(TEST_DIR)/bicg.bc.opt 2> $(RESULT_DIR)/bicg_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log

atax: $(TEST_DIR)/atax.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/atax.c -o $(TEST_DIR)/atax.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/atax.bc> $(TEST_DIR)/atax.bc.opt 2> $(RESULT_DIR)/atax_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/atax_staticSampling.cpp -O3 -o $(RESULT_DIR)/atax_staticSampling

atax_time_codeGen: $(TEST_DIR)/atax.c
	echo "atax" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/atax.bc> $(TEST_DIR)/atax.bc.opt 2> $(RESULT_DIR)/atax_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log

3mm: $(TEST_DIR)/3mm.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/3mm.c -o $(TEST_DIR)/3mm.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/3mm.bc> $(TEST_DIR)/3mm.bc.opt 2> $(RESULT_DIR)/3mm_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/3mm_staticSampling.cpp -O3 -o $(RESULT_DIR)/3mm_staticSampling

3mm_time_codeGen: $(TEST_DIR)/3mm.c
	echo "3mm" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/3mm.bc> $(TEST_DIR)/3mm.bc.opt 2> $(RESULT_DIR)/3mm_staticSampling.cpp ; } 2>> $(RESULT_DIR)/time.log

covariance: $(TEST_DIR)/covariance.c
	$(CC) -O0 -c -emit-llvm $(TEST_DIR)/covariance.c -o $(TEST_DIR)/covariance.bc
	$(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/covariance.bc> $(TEST_DIR)/covariance.bc.opt 2> $(RESULT_DIR)/covariance_staticSampling.cpp
	g++ -std=c++11 $(RESULT_DIR)/covariance_staticSampling -O3 -o $(RESULT_DIR)/covariance_staticSampling

covariance_time_codeGen: $(TEST_DIR)/covariance.c
	echo "covariance" >> $(RESULT_DIR)/time.log
	{ time $(OPT) -load /Users/dongchen/tools/llvm_xcode/Debug/lib/LLVMSymFP.dylib -symFP <$(TEST_DIR)/covariance.bc> $(TEST_DIR)/covariance.bc.opt 2> $(RESULT_DIR)/covariance_staticSampling.cpp ;} 2>> $(RESULT_DIR)/time.log

run_ss:
	$(RESULT_DIR)/stencil_staticSampling > $(RESULT_DIR)/stencil_staticSampling_result.txt &
#	$(RESULT_DIR)/correlation_staticSampling > $(RESULT_DIR)/correlation_staticSampling_result.txt
	$(RESULT_DIR)/convolution_3d_staticSampling > $(RESULT_DIR)/convolution_3d_staticSampling_result.txt &
	$(RESULT_DIR)/convolution_2d_staticSampling > $(RESULT_DIR)/convolution_2d_staticSampling_result.txt &
	$(RESULT_DIR)/syrk_staticSampling > $(RESULT_DIR)/syrk_staticSampling_result.txt &
	$(RESULT_DIR)/mvt_staticSampling > $(RESULT_DIR)/mvt_staticSampling_result.txt &
	$(RESULT_DIR)/gesummv_staticSampling > $(RESULT_DIR)/gesummv_staticSampling_result.txt &
	$(RESULT_DIR)/gemver_staticSampling > $(RESULT_DIR)/gemver_staticSampling_result.txt &
	$(RESULT_DIR)/gemm_staticSampling > $(RESULT_DIR)/gemm_staticSampling_result.txt &
	$(RESULT_DIR)/bicg_staticSampling > $(RESULT_DIR)/bicg_staticSampling_result.txt &
	$(RESULT_DIR)/atax_staticSampling > $(RESULT_DIR)/atax_staticSampling_result.txt &
	$(RESULT_DIR)/3mm_staticSampling > $(RESULT_DIR)/3mm_staticSampling_result.txt &
#	$(RESULT_DIR)/covariance_staticSampling > $(RESULT_DIR)/covariance_staticSampling_result.txt &

time_ss:
	echo "stencil" > $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/stencil_staticSampling > $(RESULT_DIR)/stencil_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log
	echo "convolution_3d" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/convolution_3d_staticSampling > $(RESULT_DIR)/convolution_3d_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log
	echo "convolution_2d" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/convolution_2d_staticSampling > $(RESULT_DIR)/convolution_2d_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log
	echo "syrk" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/syrk_staticSampling > $(RESULT_DIR)/syrk_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log 
	echo "mvt" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/mvt_staticSampling > $(RESULT_DIR)/mvt_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log
	echo "gesummv" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/gesummv_staticSampling > $(RESULT_DIR)/gesummv_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log
	echo "gemver" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/gemver_staticSampling > $(RESULT_DIR)/gemver_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log
	echo "gemm" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/gemm_staticSampling > $(RESULT_DIR)/gemm_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log
	echo "bicg" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/bicg_staticSampling > $(RESULT_DIR)/bicg_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log
	echo "atax" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/atax_staticSampling > $(RESULT_DIR)/atax_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log
	echo "3mm" >> $(RESULT_DIR)/time_ss.log
	{ time $(RESULT_DIR)/3mm_staticSampling > $(RESULT_DIR)/3mm_staticSampling_result.txt ; } 2>> $(RESULT_DIR)/time_ss.log

clean:
	rm -f polyBench/*.bc polyBench/*.opt polyBench/*.o

clean_result:
	rm polyBench/result/*.txt
