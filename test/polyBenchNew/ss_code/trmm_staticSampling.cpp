
 /* Start to analysis array index
Array index info
B.addr ((k * 1024) + j)
A.addr ((k * 1024) + i)
B.addr ((i * 1024) + j)
B.addr ((i * 1024) + j)
B.addr ((i * 1024) + j)
B.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %B
double %alpha

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------k
------Loop Bound: ((i + 1), 1024)
--------array access A.addr ((k * 1024) + i)
--------array access B.addr ((k * 1024) + j)
--------array access B.addr ((i * 1024) + j)
--------array access B.addr ((i * 1024) + j)
------array access B.addr ((i * 1024) + j)
------array access B.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
trmm */ 
