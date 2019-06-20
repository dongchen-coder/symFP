
llvm_path=/Users/dchen/tools/llvm-8.0.0/build/bin
sps_path=/Users/dchen/tools/symFP

# gen_time
'''
rm -r gen_time
mkdir gen_time

for name in 2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2d syrk trisolv trmm
do
#	time --output="./gen_time/""$name""_gen_time.log" -p sh -c "$llvm_path""opt -load ../../build/sps/libLLVMsps.so -sps  <../../test_facility/ss_bc/$name.bc> ../../test_facility/ss_bc/$name.bc.opt 2> ../../test_facility/ss_code/${name}_staticSampling.cpp"
	{ time $llvm_path/opt -load $sps_path/build/sps/libLLVMsps.so -sps  <$sps_path/test_facility/ss_bc/$name.bc> $sps_path/test_facility/ss_bc/$name.bc.opt 2> $sps_path/test_facility/ss_code/${name}_staticSampling.cpp ; } &> ./gen_time/${name}_gen_time.log
done

for i in 1 2 3 4
do
	for name in 2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2d syrk trisolv trmm
	do
		#time --output="./gen_time/""$name""_gen_time.log" -a -p sh -c "opt -load ../../build/sps/LLVMsps.so -sps  <../../test_facility/bc/$name.bc> ../../test_facility/bc/$name.bc.opt 2> ../../test_facility/ss_code/${name}_staticSampling.cpp"
		{ time $llvm_path/opt -load $sps_path/build/sps/libLLVMsps.so -sps  <$sps_path/test_facility/ss_bc/$name.bc> $sps_path/test_facility/ss_bc/$name.bc.opt 2> $sps_path/test_facility/ss_code/${name}_staticSampling.cpp ; } >> ./gen_time/${name}_gen_time.log 2>&1
	done
done
'''
# run time
#for rate in 0.1 0.02 0.01 0.005 0.002 0.001

for rate in 0.01
do

	cd $sps_path/test_run/overhead

	rm -r run_time_$rate
	mkdir run_time_$rate
	
    cd $sps_path/test_facility
    make all_gen SPSRATE=$rate

	echo "start to time"
	for name in 2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2d syrk trisolv trmm
	do
		#time --output="/sps_pldi18_aec/test_run/overhead/run_time_${rate}/${name}_run_time.log" -p sh -c "/sps_pldi18_aec/test_facility/ss_bin/${name}_staticSampling > /sps_pldi18_aec/test_facility/ss_result/${name}_staticSampling_result.txt"
		{ time $sps_path/test_facility/bin/${name}_staticSampling > $sps_path/test_facility/ss_result/${name}_staticSampling_result.txt ; } &> $sps_path/test_run/overhead_test/run_time_${rate}/${name}_run_time.log
	done
	for i in 1 2 3 4
	do
		for name in 2mm 3mm adi atax bicg cholesky correlation covariance deriche doitgen durbin fdtd_2d floyd_warshall gemm gemver gesummv gramschmidt heat_3d jacobi_1d jacobi_2d lu ludcmp mvt nussinov seidel_2d symm syr2d syrk trisolv trmm
		do
			#time --output="/sps_pldi18_aec/test_run/overhead/run_time_${rate}/${name}_run_time.log" -a -p sh -c "/sps_pldi18_aec/test_facility/ss_bin/${name}_staticSampling > /sps_pldi18_aec/test_facility/ss_result/${name}_staticSampling_result.txt"
			{ time $sps_path/test_facility/bin/${name}_staticSampling > $sps_path/test_facility/ss_result/${name}_staticSampling_result.txt ; } >> $sps_path/test_run/overhead_test/run_time_${rate}/${name}_run_time.log 2>&1
		done
	done
done



