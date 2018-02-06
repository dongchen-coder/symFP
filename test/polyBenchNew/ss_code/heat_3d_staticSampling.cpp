
 /* Start to analysis array index
Array index info
A.addr (((((i + 1) * 1024) * 1024) + (j * 1024)) + k)
A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
A.addr (((((i - 1) * 1024) * 1024) + (j * 1024)) + k)
A.addr ((((i * 1024) * 1024) + ((j + 1) * 1024)) + k)
A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
A.addr ((((i * 1024) * 1024) + ((j - 1) * 1024)) + k)
A.addr (((((i * 1024) * 1024) + (j * 1024)) + k) + 1)
A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
A.addr (((((i * 1024) * 1024) + (j * 1024)) + k) - 1)
A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
B.addr (((((i + 1) * 1024) * 1024) + (j * 1024)) + k)
B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
B.addr (((((i - 1) * 1024) * 1024) + (j * 1024)) + k)
B.addr ((((i * 1024) * 1024) + ((j + 1) * 1024)) + k)
B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
B.addr ((((i * 1024) * 1024) + ((j - 1) * 1024)) + k)
B.addr (((((i * 1024) * 1024) + (j * 1024)) + k) + 1)
B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
B.addr (((((i * 1024) * 1024) + (j * 1024)) + k) - 1)
B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %B
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (1, 10)
----i
----Loop Bound: (1, 1023)
------j
------Loop Bound: (1, 1023)
--------k
--------Loop Bound: (1, 1023)
----------array access A.addr (((((i + 1) * 1024) * 1024) + (j * 1024)) + k)
----------array access A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
----------array access A.addr (((((i - 1) * 1024) * 1024) + (j * 1024)) + k)
----------array access A.addr ((((i * 1024) * 1024) + ((j + 1) * 1024)) + k)
----------array access A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
----------array access A.addr ((((i * 1024) * 1024) + ((j - 1) * 1024)) + k)
----------array access A.addr (((((i * 1024) * 1024) + (j * 1024)) + k) + 1)
----------array access A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
----------array access A.addr (((((i * 1024) * 1024) + (j * 1024)) + k) - 1)
----------array access A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
----------array access B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
----i
----Loop Bound: (1, 1023)
------j
------Loop Bound: (1, 1023)
--------k
--------Loop Bound: (1, 1023)
----------array access B.addr (((((i + 1) * 1024) * 1024) + (j * 1024)) + k)
----------array access B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
----------array access B.addr (((((i - 1) * 1024) * 1024) + (j * 1024)) + k)
----------array access B.addr ((((i * 1024) * 1024) + ((j + 1) * 1024)) + k)
----------array access B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
----------array access B.addr ((((i * 1024) * 1024) + ((j - 1) * 1024)) + k)
----------array access B.addr (((((i * 1024) * 1024) + (j * 1024)) + k) + 1)
----------array access B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
----------array access B.addr (((((i * 1024) * 1024) + (j * 1024)) + k) - 1)
----------array access B.addr ((((i * 1024) * 1024) + (j * 1024)) + k)
----------array access A.addr ((((i * 1024) * 1024) + (j * 1024)) + k)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
heat_3d */ 
