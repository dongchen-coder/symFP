
 /* Start to analysis array index
Array index info
imgIn.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
imgIn.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgIn.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %y1
double* %imgIn
double* %y2
double* %imgOut
double %alpha

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access imgIn.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
------array access imgIn.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (1023, 0)
------array access y2.addr ((i * 1024) + j)
------array access imgIn.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access y1.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)
--j
--Loop Bound: (0, 1024)
----i
----Loop Bound: (0, 1024)
------array access imgOut.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
--j
--Loop Bound: (0, 1024)
----i
----Loop Bound: (1023, 0)
------array access y2.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access y1.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
deriche */ 
