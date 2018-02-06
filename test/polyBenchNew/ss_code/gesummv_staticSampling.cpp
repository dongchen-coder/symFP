
 /* Start to analysis array index
Array index info
y.addr i
tmp.addr i
A.addr ((i * 1024) + j)
x.addr j
tmp.addr i
tmp.addr i
B.addr ((i * 1024) + j)
x.addr j
y.addr i
y.addr i
tmp.addr i
y.addr i
y.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %n
double %alpha
double %beta
double* %A
double* %B
double* %tmp
double* %x
double* %y

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----array access tmp.addr i
----array access y.addr i
----j
----Loop Bound: (0, 1024)
------array access A.addr ((i * 1024) + j)
------array access x.addr j
------array access tmp.addr i
------array access tmp.addr i
------array access B.addr ((i * 1024) + j)
------array access x.addr j
------array access y.addr i
------array access y.addr i
----array access tmp.addr i
----array access y.addr i
----array access y.addr i

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
gesummv */ 
