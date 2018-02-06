
 /* Start to analysis array index
Array index info
A.addr (i + 1)
A.addr (i - 1)
A.addr i
B.addr i
B.addr (i - 1)
B.addr i
B.addr (i + 1)
A.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %B

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (0, 10)
----i
----Loop Bound: (1, 1023)
------array access A.addr (i - 1)
------array access A.addr i
------array access A.addr (i + 1)
------array access B.addr i
----i
----Loop Bound: (1, 1023)
------array access B.addr (i - 1)
------array access B.addr i
------array access B.addr (i + 1)
------array access A.addr i

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
jacobi_1d */ 
