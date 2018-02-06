
 /* Start to analysis array index
Array index info
A.addr ((i * 1024) + j)
u1.addr i
v1.addr j
u2.addr i
v2.addr j
A.addr ((i * 1024) + j)
x.addr i
A.addr ((j * 1024) + i)
y.addr j
x.addr i
x.addr i
z.addr i
x.addr i
w.addr i
A.addr ((i * 1024) + j)
x.addr j
w.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %n
double %alpha
double %beta
double* %A
double* %u1
double* %v1
double* %u2
double* %v2
double* %w
double* %x
double* %y
double* %z

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access A.addr ((i * 1024) + j)
------array access u1.addr i
------array access v1.addr j
------array access u2.addr i
------array access v2.addr j
------array access A.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access x.addr i
------array access A.addr ((j * 1024) + i)
------array access y.addr j
------array access x.addr i
--i
--Loop Bound: (0, 1024)
----array access x.addr i
----array access z.addr i
----array access x.addr i
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access w.addr i
------array access A.addr ((i * 1024) + j)
------array access x.addr j
------array access w.addr i

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
gemver */ 
