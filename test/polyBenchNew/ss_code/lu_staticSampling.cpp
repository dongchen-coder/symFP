
 /* Start to analysis array index
Array index info
A.addr ((k * 1024) + j)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((j * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
A.addr ((k * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, i)
------k
------Loop Bound: (0, j)
--------array access A.addr ((i * 1024) + k)
--------array access A.addr ((k * 1024) + j)
--------array access A.addr ((i * 1024) + j)
--------array access A.addr ((i * 1024) + j)
------array access A.addr ((j * 1024) + j)
------array access A.addr ((i * 1024) + j)
------array access A.addr ((i * 1024) + j)
----j
----Loop Bound: (i, 1024)
------k
------Loop Bound: (0, i)
--------array access A.addr ((i * 1024) + k)
--------array access A.addr ((k * 1024) + j)
--------array access A.addr ((i * 1024) + j)
--------array access A.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
lu */ 
