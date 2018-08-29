rm -r ss_result_*
for i in 0.05 0.02 0.01 0.005 0.002 0.001 
do
	
	cd ../../test_facility
	make all_gen SPSRATE=$i
	make all_run
	cd ../test_run/precision_test
	cp -r ../../test_facility/ss_result "ss_result_$i"
done
