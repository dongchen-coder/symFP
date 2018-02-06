
 /* Start to analysis array index
Array index info
A.addr (((i * 1024) + j) - 1)
A.addr ((i * 1024) + j)
A.addr (((i * 1024) + 1) + j)
A.addr (((1 + i) * 1024) + j)
A.addr (((i - 1) * 1024) + j)
B.addr ((i * 1024) + j)
B.addr ((i * 1024) + j)
B.addr (((i * 1024) + j) - 1)
B.addr (((i * 1024) + 1) + j)
B.addr (((1 + i) * 1024) + j)
B.addr (((i - 1) * 1024) + j)
A.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %B

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (0, 10)
----i
----Loop Bound: (1, 1023)
------j
------Loop Bound: (1, 1023)
--------array access A.addr ((i * 1024) + j)
--------array access A.addr (((i * 1024) + j) - 1)
--------array access A.addr (((i * 1024) + 1) + j)
--------array access A.addr (((1 + i) * 1024) + j)
--------array access A.addr (((i - 1) * 1024) + j)
--------array access B.addr ((i * 1024) + j)
----i
----Loop Bound: (1, 1023)
------j
------Loop Bound: (1, 1023)
--------array access B.addr ((i * 1024) + j)
--------array access B.addr (((i * 1024) + j) - 1)
--------array access B.addr (((i * 1024) + 1) + j)
--------array access B.addr (((1 + i) * 1024) + j)
--------array access B.addr (((i - 1) * 1024) + j)
--------array access A.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
jacobi_2d */ 
