
 /* Start to analysis array index
Array index info
C.addr ((i * 256) + j)
C.addr ((i * 256) + j)
A.addr ((i * 256) + k)
B.addr ((k * 256) + j)
C.addr ((i * 256) + j)
C.addr ((i * 256) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %ni
i32 %nj
i32 %nk
double %alpha
double %beta
double* %A
double* %B
double* %C

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 256)
----j
----Loop Bound: (0, 256)
------array access C.addr ((i * 256) + j)
------array access C.addr ((i * 256) + j)
------k
------Loop Bound: (0, 256)
--------array access A.addr ((i * 256) + k)
--------array access B.addr ((k * 256) + j)
--------array access C.addr ((i * 256) + j)
--------array access C.addr ((i * 256) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
gemm */ 
