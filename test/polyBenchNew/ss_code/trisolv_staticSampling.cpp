
 /* Start to analysis array index
Array index info
b.addr i
L.addr ((i * 1024) + j)
x.addr j
x.addr i
x.addr i
x.addr i
x.addr i
L.addr ((i * 1024) + i)
x.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %x
double* %b
double* %L

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----array access b.addr i
----array access x.addr i
----j
----Loop Bound: (0, i)
------array access L.addr ((i * 1024) + j)
------array access x.addr j
------array access x.addr i
------array access x.addr i
----array access x.addr i
----array access L.addr ((i * 1024) + i)
----array access x.addr i

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
trisolv */ 
