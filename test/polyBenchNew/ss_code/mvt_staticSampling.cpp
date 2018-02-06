
 /* Start to analysis array index
Array index info
x1.addr i
a.addr ((i * 1024) + j)
y1.addr j
x1.addr i
x2.addr i
a.addr ((j * 1024) + i)
y2.addr j
x2.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %n
double* %a
double* %x1
double* %x2
double* %y1
double* %y2

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access x1.addr i
------array access a.addr ((i * 1024) + j)
------array access y1.addr j
------array access x1.addr i
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access x2.addr i
------array access a.addr ((j * 1024) + i)
------array access y2.addr j
------array access x2.addr i

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
mvt */ 
