#===- runtime/Makefile -------------------------*- Makefile -*-===##

# Relative path to the top of the source tree.
#
LEVEL=..

#CLANG_BIN_DIR=/Users/dchen/tools/llvm-7.0.1/build/Debug/bin
#LLVM_BIN_DIR=/Users/dchen/tools/llvm-7.0.1/build/Debug/bin
#TEST_DIR=/Users/dchen/tools/symFP/test_facility
#LIB_DIR=/Users/dchen/tools/llvm-7.0.1/build/Debug/lib

LLVM_BIN_DIR=/localdisk/fliu14/llvm/build/bin
CLANG_BIN_DIR=/localdisk/fliu14/llvm/build/bin
PAPI_LIB_DIR=/localdisk/local/include
PAPI_LIB=/localdisk/local/lib/libpapi.a
TEST_DIR=.
#LIB_DIR=/Users/dchen/tools/llvm-8.0.0/build/lib
LIB_DIR=../build/lib

SPSRATE?=0.1

SRC_DIR=$(TEST_DIR)/ss_src
BC_DIR=$(TEST_DIR)/ss_bc
IR_DIR=$(TEST_DIR)/ss_ir
RESULT_DIR=$(TEST_DIR)/ss_result
SS_DIR=$(TEST_DIR)/ss_code
BIN_DIR=$(TEST_DIR)/bin

TRACE_SRC_DIR=$(TEST_DIR)/trace_src
TRACE_RESULT_DIR=$(TEST_DIR)/trace_result

TRACE_REF_SRC_DIR=$(TEST_DIR)/ref_trace_src
TRACE_REF_RESULT_DIR=$(TEST_DIR)/ref_trace_result

TRACE_REF_OSL_SRC_DIR=$(TEST_DIR)/ref_trace_osl_src
TRACE_REF_OSL_RESULT_DIR=$(TEST_DIR)/ref_trace_osl_result

OPT= $(LLVM_BIN_DIR)/opt
LLC= $(LLVM_BIN_DIR)/llc
CC= $(CLANG_BIN_DIR)/clang
DIS= $(LLVM_BIN_DIR)/llvm-dis

CHK_SIZE=4

all_bench = stencil 2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2k syrk trisolv trmm

constant_bound_bench = mm 3mm mvt atax bicg gemm gemver gesummv convolution2d convolution3d syrk syr2k stencil doitgen # deriche doitgen fdtd_2d floyd_warshall heat_3d jacobi_1d jacobi_2d seidel_2d adi

sym_bound_bench = cholesky correlation covariance durbin lu ludcmp nussinov symm trisolv trmm

cl_bench = 2mm_kernel 3mm_kernel adi_kernel atax_kernel bicg_kernel convolution_2d_kernel convolution_3d_kernel correlation_kernel covariance_kernel doitgen_kernel fdtd_kernel gemm_kernel gemver_kernel gesummv_kernel gramschmidt_kernel jacobi_1d_kernel jacobi_2d_kernel lu_kernel mvt_kernel syr2k_kernel syrk_kernel

test_bench = test1 test2 # test3 test4 test5

check_bench = 2mm 3mm mvt atax bicg gemm gemver gesummv convolution2d convolution3d syrk syr2k stencil# 
#check_bench = stencil convolution2d convolution3d
#check_bench = atax
#check_bench = gesummv
#check_bench = heat_3d
#check_bench = cholesky
#check_bench = correlation
#check_bench = covariance
#check_bench = durbin
#check_bench = lu
#check_bench = ludcmp
#check_bench = nussinov
#check_bench = symm
#check_bench = syr2k
#check_bench = syrk
#check_bench = trisolv
#check_bench = trmm
#check_bench = 3mm
#check_bench = mm

all_gen:
	$(foreach name, $(all_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(all_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(all_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/sequential/$(name)_staticSampling.cpp ;) 
	$(foreach name, $(all_bench), g++ -std=c++11 -g $(SS_DIR)/sequential/$(name)_staticSampling.cpp -O3 -pthread -DPAPI_TIMER  -o $(BIN_DIR)/sequential/$(name)_staticSampling ;)

parallel_gen:
	$(foreach name, $(all_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(all_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(all_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/parallel/$(name)_staticSampling.cpp ;)
	$(foreach name, $(all_bench), g++ -std=c++11 $(SS_DIR)/parallel/$(name)_staticSampling.cpp -O3 -DTHREAD_NUM=8 -DBIN_SIZE=64 -o $(BIN_DIR)/parallel/$(name)_staticSampling ;)

all_run:
	$(foreach name, $(all_bench), $(BIN_DIR)/$(name)_staticSampling > $(RESULT_DIR)/$(name)_staticSampling_result.txt ;)

constant_gen:
	$(foreach name, $(constant_bound_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(constant_bound_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(constant_bound_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/sequential/$(name)_staticSampling.cpp ;)
	$(foreach name, $(constant_bound_bench), g++ -std=c++11 $(SS_DIR)/sequential/$(name)_staticSampling.cpp -O3 -o $(BIN_DIR)/sequential/$(name)_staticSampling ;)

constant_parallel_gen:
	$(foreach name, $(constant_bound_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -DSMALL_DATASET -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(constant_bound_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(constant_bound_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/parallel/$(name)_staticSampling.cpp ;)
	$(foreach name, $(constant_bound_bench), g++ -std=c++11 $(SS_DIR)/parallel/$(name)_staticSampling.cpp -O3 -o $(BIN_DIR)/parallel/$(name)_staticSampling ;)

debug:
	$(foreach name, $(check_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -DDEBUG -DMINI_DATASET -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(check_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=1 <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/parallel/$(name)_staticSampling.cpp ;)
	$(foreach name, $(check_bench), g++ -std=c++11 $(SS_DIR)/parallel/$(name)_staticSampling.cpp -O3 -DTHREAD_NUM=2 -DCHUNK_SIZE=2 -DDEBUG -o $(BIN_DIR)/parallel/$(name)_staticSampling ;)

constant_run:
	$(foreach name, $(constant_bound_bench), $(BIN_DIR)/$(name)_staticSampling > $(RESULT_DIR)/$(name)_staticSampling_result.txt &)

sym_gen:
	$(foreach name, $(sym_bound_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(sym_bound_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(sym_bound_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)
	$(foreach name, $(sym_bound_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)

sym_run:
	$(foreach name, $(sym_bound_bench), $(BIN_DIR)/$(name)_staticSampling > $(RESULT_DIR)/$(name)_staticSampling_result.txt &)

cl_gen:
	$(foreach name, $(cl_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).cl -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(cl_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(cl_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)
	$(foreach name, $(cl_bench), g++ -std=c++11 $(SS_DIR)/$(name)_staticSampling.cpp -O3 -o $(BIN_DIR)/$(name)_staticSampling ;)

check_gen:
	$(foreach name, $(check_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -DSMALL_DATASET -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(check_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(check_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)
	$(foreach name, $(check_bench), g++ -g -std=c++11 -lpthread $(SS_DIR)/$(name)_staticSampling.cpp -O0 -DBIN_SIZE=64 -DTHREAD_NUM=8 -DCHUNK_SIZE=$(CHK_SIZE) -o $(BIN_DIR)/$(name)_staticSampling ;)

check_run:
	$(foreach name, $(check_bench), $(BIN_DIR)/$(name)_staticSampling > $(RESULT_DIR)/$(name)_staticSampling_result.txt ;)

trace_run:
	$(foreach name, $(check_bench), g++ -std=c++11 $(TRACE_SRC_DIR)/$(name).cpp -O3 -DMINI_DATASET -o $(BIN_DIR)/$(name)_trace ;)
	$(foreach name, $(check_bench), $(BIN_DIR)/$(name)_trace > $(TRACE_RESULT_DIR)/$(name)_trace_result.txt ;)

trace_ref_run:
	$(foreach name, $(all_bench), g++ -std=c++11 $(TRACE_REF_SRC_DIR)/$(name)_ref_arr.cpp -O3 -o $(BIN_DIR)/$(name)_ref_arr_trace ;)
	$(foreach name, $(all_bench), $(BIN_DIR)/$(name)_ref_arr_trace > $(TRACE_REF_RESULT_DIR)/$(name)_ref_arr_trace_result.txt ;)

trace_ref_osl_run:
	$(foreach name, $(all_bench), g++ -std=c++11 $(TRACE_REF_OSL_SRC_DIR)/$(name)_ref_arr.cpp -O3 -I$(PAPI_LIB_DIR) $(PAPI_LIB) -DPAPI_TIMER -o $(BIN_DIR)/$(name)_ref_osl_trace ;)
	# $(foreach name, $(all_bench), $(BIN_DIR)/$(name)_ref_osl_trace > $(TRACE_REF_OSL_RESULT_DIR)/$(name)_ref_osl_trace_result.txt ;)

test_trace_run:
	$(foreach name, $(test_bench), g++ -std=c++11 $(TRACE_SRC_DIR)/$(name).cpp -O3 -DMINI_DATASET -o $(BIN_DIR)/$(name)_trace ;)

test_parallel_gen:
	$(foreach name, $(test_bench), $(CC) -O0 -c -emit-llvm $(SRC_DIR)/$(name).c -DMINI_DATASET -o $(BC_DIR)/$(name).bc ;)
	$(foreach name, $(test_bench), $(DIS) $(BC_DIR)/$(name).bc -o $(IR_DIR)/$(name).ir ;)
	$(foreach name, $(test_bench), $(OPT) -load $(LIB_DIR)/libLLVMSPS.so -sps -spsrate=$(SPSRATE) <$(BC_DIR)/$(name).bc> $(BC_DIR)/$(name).bc.opt 2> $(SS_DIR)/$(name)_staticSampling.cpp ;)
	$(foreach name, $(test_bench), g++ -g -std=c++11 $(SS_DIR)/$(name)_staticSampling.cpp -O3 -DBIN_SIZE=64 -DTHREAD_NUM=8 -DCHUNK_SIZE=$(CHK_SIZE) -o $(BIN_DIR)/$(name)_staticSampling ;)

clean:	
	rm $(BIN_DIR)/sequential/*_staticSampling
	rm $(BIN_DIR)/parallel/*_staticSampling
	rm $(BIN_DIR)/*_staticSampling

veryclean:
	rm $(BC_DIR)/*.bc
	rm $(IR_DIR)/*.ir
	rm $(SS_DIR)/sequential/*_staticSampling.cpp
	rm $(SS_DIR)/parallel/*_staticSampling.cpp
	rm $(BIN_DIR)/sequential/*_staticSampling
	rm $(BIN_DIR)/parallel/*_staticSampling
	rm $(BIN_DIR)/*_staticSampling
