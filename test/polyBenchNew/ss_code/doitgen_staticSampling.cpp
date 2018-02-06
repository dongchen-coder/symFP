
 /* Start to analysis array index
Array index info
A.addr ((((r * 256) * 256) + (q * 256)) + s)
sum.addr p
C4.addr ((s * 256) + p)
sum.addr p
sum.addr p
sum.addr p
A.addr ((((r * 256) * 256) + (q * 256)) + p)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %sum
double* %A
double* %C4

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--r
--Loop Bound: (0, 256)
----q
----Loop Bound: (0, 256)
------p
------Loop Bound: (0, 256)
--------array access sum.addr p
--------s
--------Loop Bound: (0, 256)
----------array access A.addr ((((r * 256) * 256) + (q * 256)) + s)
----------array access C4.addr ((s * 256) + p)
----------array access sum.addr p
----------array access sum.addr p
------p
------Loop Bound: (0, 256)
--------array access sum.addr p
--------array access A.addr ((((r * 256) * 256) + (q * 256)) + p)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
doitgen */ 
