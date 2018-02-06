
 /* Start to analysis array index
Array index info
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
A.addr ((k * 1024) + j)
A.addr ((j * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
A.addr ((k * 1024) + j)
A.addr ((i * 1024) + j)
b.addr i
A.addr ((i * 1024) + j)
y.addr j
y.addr i
y.addr i
A.addr ((i * 1024) + j)
x.addr j
A.addr ((i * 1024) + i)
x.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %b
double* %y
double* %x

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, i)
------array access A.addr ((i * 1024) + j)
------k
------Loop Bound: (0, j)
--------array access A.addr ((i * 1024) + k)
--------array access A.addr ((k * 1024) + j)
------array access A.addr ((j * 1024) + j)
------array access A.addr ((i * 1024) + j)
----j
----Loop Bound: (i, 1024)
------array access A.addr ((i * 1024) + j)
------k
------Loop Bound: (0, i)
--------array access A.addr ((i * 1024) + k)
--------array access A.addr ((k * 1024) + j)
------array access A.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
----array access b.addr i
----j
----Loop Bound: (0, i)
------array access A.addr ((i * 1024) + j)
------array access y.addr j
----array access y.addr i
--i
--Loop Bound: (1023, 0)
----array access y.addr i
----j
----Loop Bound: ((i + 1), 1024)
------array access A.addr ((i * 1024) + j)
------array access x.addr j
----array access A.addr ((i * 1024) + i)
----array access x.addr i

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
ludcmp */ 
