
 /* Start to analysis array index
Array index info
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + k)
R.addr ((k * 1024) + k)
A.addr ((i * 1024) + k)
R.addr ((k * 1024) + k)
Q.addr ((i * 1024) + k)
R.addr ((k * 1024) + j)
Q.addr ((i * 1024) + k)
A.addr ((i * 1024) + j)
R.addr ((k * 1024) + j)
R.addr ((k * 1024) + j)
A.addr ((i * 1024) + j)
Q.addr ((i * 1024) + k)
R.addr ((k * 1024) + j)
A.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %R
double* %Q

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--k
--Loop Bound: (0, 1024)
----i
----Loop Bound: (0, 1024)
------array access A.addr ((i * 1024) + k)
------array access A.addr ((i * 1024) + k)
----array access R.addr ((k * 1024) + k)
----i
----Loop Bound: (0, 1024)
------array access A.addr ((i * 1024) + k)
------array access R.addr ((k * 1024) + k)
------array access Q.addr ((i * 1024) + k)
----j
----Loop Bound: ((k + 1), 1024)
------array access R.addr ((k * 1024) + j)
------i
------Loop Bound: (0, 1024)
--------array access Q.addr ((i * 1024) + k)
--------array access A.addr ((i * 1024) + j)
--------array access R.addr ((k * 1024) + j)
--------array access R.addr ((k * 1024) + j)
------i
------Loop Bound: (0, 1024)
--------array access A.addr ((i * 1024) + j)
--------array access Q.addr ((i * 1024) + k)
--------array access R.addr ((k * 1024) + j)
--------array access A.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
gramschmidt */ 
