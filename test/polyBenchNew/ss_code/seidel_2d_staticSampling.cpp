
 /* Start to analysis array index
Array index info
A.addr (((i - 1) * 1024) + j)
A.addr ((((i - 1) * 1024) + j) - 1)
A.addr ((((i - 1) * 1024) + j) + 1)
A.addr (((i * 1024) + j) - 1)
A.addr ((i * 1024) + j)
A.addr (((i * 1024) + j) + 1)
A.addr ((((i + 1) * 1024) + j) - 1)
A.addr (((i + 1) * 1024) + j)
A.addr ((((i + 1) * 1024) + j) + 1)
A.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (0, 9)
----i
----Loop Bound: (1, 1022)
------j
------Loop Bound: (1, 1022)
--------array access A.addr ((((i - 1) * 1024) + j) - 1)
--------array access A.addr (((i - 1) * 1024) + j)
--------array access A.addr ((((i - 1) * 1024) + j) + 1)
--------array access A.addr (((i * 1024) + j) - 1)
--------array access A.addr ((i * 1024) + j)
--------array access A.addr (((i * 1024) + j) + 1)
--------array access A.addr ((((i + 1) * 1024) + j) - 1)
--------array access A.addr (((i + 1) * 1024) + j)
--------array access A.addr ((((i + 1) * 1024) + j) + 1)
--------array access A.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
seidel_2d */ 
