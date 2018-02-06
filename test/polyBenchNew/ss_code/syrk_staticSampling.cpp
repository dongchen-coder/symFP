
 /* Start to analysis array index
Array index info
C.addr ((i * 256) + j)
C.addr ((i * 256) + j)
A.addr ((i * 256) + k)
A.addr ((j * 256) + k)
C.addr ((i * 256) + j)
C.addr ((i * 256) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %ni
i32 %nj
double %alpha
double %beta
double* %A
double* %C

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 256)
----j
----Loop Bound: (0, i)
------array access C.addr ((i * 256) + j)
------array access C.addr ((i * 256) + j)
------k
------Loop Bound: (0, 256)
--------j
--------Loop Bound: (0, i)
----------array access A.addr ((i * 256) + k)
----------array access A.addr ((j * 256) + k)
----------array access C.addr ((i * 256) + j)
----------array access C.addr ((i * 256) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
syrk */ 
