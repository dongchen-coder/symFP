
 /* Start to analysis array index
Array index info
r.addr 0
y.addr 0
r.addr 0
r.addr ((k - i) - 1)
y.addr i
r.addr k
y.addr i
y.addr ((k - i) - 1)
z.addr i
z.addr i
y.addr i
y.addr k

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %y
double* %r
double* %z

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--array access r.addr 0
--array access y.addr 0
--array access r.addr 0
--k
--Loop Bound: (1, 1024)
----i
----Loop Bound: (0, k)
------array access r.addr ((k - i) - 1)
------array access y.addr i
----array access r.addr k
----i
----Loop Bound: (0, k)
------array access y.addr i
------array access y.addr ((k - i) - 1)
------array access z.addr i
----i
----Loop Bound: (0, k)
------array access z.addr i
------array access y.addr i
----array access y.addr k

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
durbin */ 
