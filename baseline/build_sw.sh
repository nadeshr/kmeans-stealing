make clean 
make
make hw-sim 
make run | tee exec_trace 
grep -B 38 "Execution time" exec_trace
