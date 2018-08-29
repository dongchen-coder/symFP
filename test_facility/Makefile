#===- runtime/Makefile -------------------------*- Makefile -*-===##

# Relative path to the top of the source tree.
#
LEVEL=..

#LLVM_BIN_DIR=/u/cs255/build-llvm-38/bin
CLANG_BIN_DIR=/Users/dongchen/tools/llvm-6.0.0.build-xcode/Debug/bin
LLVM_BIN_DIR=/Users/dongchen/tools/llvm-6.0.0.build-xcode/Debug/bin
TEST_DIR=/Users/dongchen/tools/llvm-6.0.0.src/lib/Transforms/StaticSampling/test_facility
LIB_DIR=/Users/dongchen/tools/llvm-6.0.0.build-xcode/Debug/lib

SPSRATE?=0.01

SRC_DIR=$(TEST_DIR)/ss_src
BC_DIR=$(TEST_DIR)/ss_bc
IR_DIR=$(TEST_DIR)/ss_ir
RESULT_DIR=$(TEST_DIR)/ss_result
SS_DIR=$(TEST_DIR)/ss_code
BIN_DIR=$(TEST_DIR)/bin

TRACE_SRC_DIR=$(TEST_DIR)/trace_src
TRACE_RESULT_DIR=$(TEST_DIR)/trace_result

OPT= $(LLVM_BIN_DIR)/opt
LLC= $(LLVM_BIN_DIR)/llc
CC= $(CLANG_BIN_DIR)/clang
DIS= $(LLVM_BIN_DIR)/llvm-dis

all_bench = 2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2d syrk trisolv trmm

constant_bound_bench = 2mm 3mm adi atax bicg deriche doitgen fdtd_2d floyd_warshall gemm gemver gesummv heat_3d jacobi_1d jacobi_2d mvt seidel_2d

sym_bound_bench = cholesky correlation covariance durbin lu ludcmp nussinov symm syr2d syrk trisolv trmm

cl_bench = 2mm_kernel 3mm_kernel adi_kernel atax_kernel bicg_kernel convolution_2d_kernel convolution_3d_kernel correlation_kernel covariance_kernel doitgen_kernel fdtd_kernel gemm_kernel gemver_kernel gesummv_kernel gramschmidt_kernel jacobi_1d_kernel jacobi_2d_kernel lu_kernel mvt_kernel syr2k_kernel syrk_kernel

check_bench = stencil
#check_bench = atax
#check_bench = 2mm
#check_bench = heat_3d
#check_bench = cholesky
#check_bench = correlation
#check_bench = covariance
#check_bench = durbin
#check_bench = lu
#check_bench = ludcmp
#check_bench = nussinov
#check_bench = symm
#check_bench = syr2d
#check_bench = syrk
#check_bench = trisolv
#check_bench = trmm


all_gen:
	$(foreach name, $(all_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(all_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(all_bench), $(OPT) -load $(LIB_DIR)/LLVMStaticSampling.dylib -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)
	$(foreach name, $(all_bench), g++ -std=c++11 $(SS_DIR)/$(name)_staticSampling.cpp -O3 -o $(BIN_DIR)/$(name)_staticSampling ;)

all_run:
	$(foreach name, $(all_bench), $(BIN_DIR)/$(name)_staticSampling > $(RESULT_DIR)/$(name)_staticSampling_resut.txt &)

constant_gen:
	$(foreach name, $(constant_bound_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(constant_bound_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(constant_bound_bench), $(OPT) -load $(LIB_DIR)/LLVMStaticSampling.dylib -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)
	$(foreach name, $(constant_bound_bench), g++ -std=c++11 $(SS_DIR)/$(name)_staticSampling.cpp -O3 -o $(BIN_DIR)/$(name)_staticSampling ;)

constant_run:
	$(foreach name, $(constant_bound_bench), $(BIN_DIR)/$(name)_staticSampling > $(RESULT_DIR)/$(name)_staticSampling_resut.txt &)

sym_gen:
	$(foreach name, $(sym_bound_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(sym_bound_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(sym_bound_bench), $(OPT) -load $(LIB_DIR)/LLVMStaticSampling.dylib -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)
	$(foreach name, $(sym_bound_bench), g++ -std=c++11 $(SS_DIR)/$(name)_staticSampling.cpp -O3 -o $(BIN_DIR)/$(name)_staticSampling ;)

sym_run:
	$(foreach name, $(sym_bound_bench), $(BIN_DIR)/$(name)_staticSampling > $(RESULT_DIR)/$(name)_staticSampling_resut.txt &)

cl_gen:
	$(foreach name, $(cl_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).cl -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(cl_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(cl_bench), $(OPT) -load $(LIB_DIR)/LLVMStaticSampling.dylib -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)
	$(foreach name, $(cl_bench), g++ -std=c++11 $(SS_DIR)/$(name)_staticSampling.cpp -O3 -o $(BIN_DIR)/$(name)_staticSampling ;)

check_gen:
	$(foreach name, $(check_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(check_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(check_bench), $(OPT) -load $(LIB_DIR)/LLVMStaticSampling.dylib -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)
	$(foreach name, $(check_bench), g++ -std=c++11 $(SS_DIR)/$(name)_staticSampling.cpp -O3 -o $(BIN_DIR)/$(name)_staticSampling ;)

check_run:
	$(foreach name, $(check_bench), $(BIN_DIR)/$(name)_staticSampling > $(RESULT_DIR)/$(name)_staticSampling_resut.txt ;)

trace_run:
	$(foreach name, $(all_bench), g++ -std=c++11 $(TRACE_SRC_DIR)/$(name).cpp -O3 -o $(BIN_DIR)/$(name)_trace ;)
	$(foreach name, $(all_bench), $(BIN_DIR)/$(name)_trace > $(TRACE_RESULT_DIR)/$(name)_trace_result.txt ;)


clean:
	rm $(BC_DIR)/*.bc
	rm $(IR_DIR)/*.ir
	rm $(SS_DIR)/*_staticSampling.cpp
	rm $(BIN_DIR)/*_staticSampling

