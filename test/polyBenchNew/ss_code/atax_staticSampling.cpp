
 /* Start to analysis array index
Array index info
tmp.addr i
y.addr i
tmp.addr i
A.addr ((i * 1024) + j)
x.addr j
tmp.addr i
y.addr j
A.addr ((i * 1024) + j)
tmp.addr i
y.addr j

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %nx
i32 %ny
double* %A
double* %x
double* %y
double* %tmp

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----array access y.addr i
--i
--Loop Bound: (0, 1024)
----array access tmp.addr i
----j
----Loop Bound: (0, 1024)
------array access tmp.addr i
------array access A.addr ((i * 1024) + j)
------array access x.addr j
------array access tmp.addr i
----j
----Loop Bound: (0, 1024)
------array access y.addr j
------array access A.addr ((i * 1024) + j)
------array access tmp.addr i
------array access y.addr j

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
atax */ 
