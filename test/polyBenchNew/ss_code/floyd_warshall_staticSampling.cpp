
 /* Start to analysis array index
Array index info
path.addr ((k * 1024) + j)
path.addr ((i * 1024) + j)
path.addr ((i * 1024) + k)
path.addr ((i * 1024) + j)
path.addr ((i * 1024) + k)
path.addr ((k * 1024) + j)
path.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %path

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--k
--Loop Bound: (0, 1024)
----i
----Loop Bound: (0, 1024)
------j
------Loop Bound: (0, 1024)
--------array access path.addr ((i * 1024) + j)
--------array access path.addr ((i * 1024) + k)
--------array access path.addr ((k * 1024) + j)
--------array access path.addr ((i * 1024) + j)
--------array access path.addr ((i * 1024) + k)
--------array access path.addr ((k * 1024) + j)
--------array access path.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
floyd_warshall */ 
