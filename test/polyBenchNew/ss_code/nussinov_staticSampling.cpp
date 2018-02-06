
 /* Start to analysis array index
Array index info
table.addr ((i * 1024) + j)
table.addr (((i * 1024) + j) - 1)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((((i + 1) * 1024) + j) - 1)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + k)
table.addr (((k + 1) * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr (((i + 1) * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((i * 1024) + j)
table.addr ((((i + 1) * 1024) + j) - 1)
seq.addr i
seq.addr j

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %table
double* %seq

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (1023, 0)
----j
----Loop Bound: ((i + 1), 1024)
------array access table.addr ((i * 1024) + j)
------array access table.addr (((i * 1024) + j) - 1)
------array access table.addr ((i * 1024) + j)
------array access table.addr ((i * 1024) + j)
------array access table.addr (((i + 1) * 1024) + j)
------array access table.addr ((i * 1024) + j)
------array access table.addr ((i * 1024) + j)
------array access table.addr ((((i + 1) * 1024) + j) - 1)
------array access seq.addr i
------array access seq.addr j
------array access table.addr ((i * 1024) + j)
------array access table.addr ((i * 1024) + j)
------array access table.addr ((((i + 1) * 1024) + j) - 1)
------array access table.addr ((i * 1024) + j)
------k
------Loop Bound: ((i + 1), j)
--------array access table.addr ((i * 1024) + j)
--------array access table.addr ((i * 1024) + k)
--------array access table.addr (((k + 1) * 1024) + j)
--------array access table.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
nussinov */ 
