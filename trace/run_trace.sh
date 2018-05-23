mkdir time_logs

rm -r ./polyBench/time_log
mkdir ./polyBench/time_log
echo "#define ORG" > ./utility/data_size.h
make
make time_trace
make time_trace
make time_trace
make time_trace
make time_trace
rm -r ./time_logs/org_log
mkdir ./time_logs/org_log
mv ./polyBench/time_log/* ./time_logs/org_log

rm -r ./polyBench/time_log
mkdir ./polyBench/time_log
echo "#define TX" > ./utility/data_size.h
make
make time_trace
make time_trace
make time_trace
make time_trace
make time_trace
rm -r ./time_logs/tx_log
mkdir ./time_logs/tx_log
mv ./polyBench/time_log/* ./time_logs/tx_log

rm -r ./polyBench/time_log
mkdir ./polyBench/time_log
echo "#define FX" > ./utility/data_size.h
make
make time_trace
make time_trace
make time_trace
make time_trace
make time_trace
rm -r ./time_logs/fx_log
mkdir ./time_logs/fx_log
mv ./polyBench/time_log/* ./time_logs/fx_log

rm -r ./polyBench/time_log
mkdir ./polyBench/time_log
echo "#define EX" > ./utility/data_size.h
make
make time_trace
make time_trace
make time_trace
make time_trace
make time_trace
rm -r ./time_logs/ex_log
mkdir ./time_logs/ex_log
mv ./polyBench/time_log/* ./time_logs/ex_log

