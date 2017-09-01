
TEST_DIR=/Users/dongchen/tools/llvm-4.0.0.src/lib/Transforms/SymFP/trace/polyBench

polyBench: stencil convolution_3d convolution_2d syrk mvt gesummv gemver gemm bicg atax 3mm

stencil:
	g++ $(TEST_DIR)/stencil.cpp -O3 -o $(TEST_DIR)/stencil_trace

convolution_3d:
	g++ $(TEST_DIR)/convolution_3d.cpp -O3 -o $(TEST_DIR)/convolution_3d_trace

convolution_2d:
	g++ $(TEST_DIR)/convolution_2d.cpp -O3 -o $(TEST_DIR)/convolution_2d_trace

syrk:
	g++ $(TEST_DIR)/syrk.cpp -O3 -o $(TEST_DIR)/syrk_trace

mvt:
	g++ $(TEST_DIR)/mvt.cpp -O3 -o $(TEST_DIR)/mvt_trace

gesummv:
	g++ $(TEST_DIR)/gesummv.cpp -O3 -o $(TEST_DIR)/gesummv_trace

gemver:
	g++ $(TEST_DIR)/gemver.cpp -O3 -o $(TEST_DIR)/gemver_trace

gemm:
	g++ $(TEST_DIR)/gemm.cpp -O3 -o $(TEST_DIR)/gemm_trace

bicg:
	g++ $(TEST_DIR)/bicg.cpp -O3 -o $(TEST_DIR)/bicg_trace

atax:
	g++ $(TEST_DIR)/atax.cpp -O3 -o $(TEST_DIR)/atax_trace

3mm:
	g++ $(TEST_DIR)/3mm.cpp -O3 -o $(TEST_DIR)/3mm_trace

run_trace:
	$(TEST_DIR)/stencil_trace > $(TEST_DIR)/result/stencil_trace_result.txt
	$(TEST_DIR)/convolution_3d_trace > $(TEST_DIR)/result/convolution_3d_trace_result.txt
	$(TEST_DIR)/convolution_2d_trace > $(TEST_DIR)/result/convolution_2d_trace_result.txt
	$(TEST_DIR)/syrk_trace > $(TEST_DIR)/result/syrk_trace_result.txt
	$(TEST_DIR)/mvt_trace > $(TEST_DIR)/result/mvt_trace_result.txt
	$(TEST_DIR)/gesummv_trace > $(TEST_DIR)/result/gesummv_trace_result.txt
	$(TEST_DIR)/gemver_trace > $(TEST_DIR)/result/gemver_trace_result.txt
	$(TEST_DIR)/gemm_trace > $(TEST_DIR)/result/gemm_trace_result.txt
	$(TEST_DIR)/bicg_trace > $(TEST_DIR)/result/bicg_trace_result.txt
	$(TEST_DIR)/atax_trace > $(TEST_DIR)/result/atax_trace_result.txt
	$(TEST_DIR)/3mm_trace > $(TEST_DIR)/result/3mm_trace_result.txt

