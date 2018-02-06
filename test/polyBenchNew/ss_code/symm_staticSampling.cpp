
 /* Start to analysis array index
Array index info
B.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
C.addr ((k * 1024) + j)
C.addr ((k * 1024) + j)
B.addr ((k * 1024) + j)
A.addr ((i * 1024) + k)
C.addr ((i * 1024) + j)
B.addr ((i * 1024) + j)
A.addr ((i * 1024) + i)
C.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %B
double* %C
double %alpha
double %beta

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------k
------Loop Bound: (0, i)
--------array access B.addr ((i * 1024) + j)
--------array access A.addr ((i * 1024) + k)
--------array access C.addr ((k * 1024) + j)
--------array access C.addr ((k * 1024) + j)
--------array access B.addr ((k * 1024) + j)
--------array access A.addr ((i * 1024) + k)
------array access C.addr ((i * 1024) + j)
------array access B.addr ((i * 1024) + j)
------array access A.addr ((i * 1024) + i)
------array access C.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
symm */ 
