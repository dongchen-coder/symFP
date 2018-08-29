cd ../../test_facility
make trace_run
cd ../test_run/precision_test
rm -r trace_result
cp -r ../../test_facility/trace_result trace_result

