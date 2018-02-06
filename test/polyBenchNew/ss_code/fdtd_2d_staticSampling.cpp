
 /* Start to analysis array index
Array index info
ex.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)
ex.addr (((i * 1024) + j) + 1)
ex.addr ((i * 1024) + j)
ey.addr (((i + 1) * 1024) + j)
ey.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)
ey.addr (0 + j)
_fict_.addr t
ey.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)
hz.addr (((i - 1) * 1024) + j)
ey.addr ((i * 1024) + j)
ex.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)
hz.addr (((i * 1024) + j) - 1)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %_fict_
double* %ey
double* %ex
double* %hz

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (0, 10)
----j
----Loop Bound: (0, 1024)
------array access _fict_.addr t
------array access ey.addr (0 + j)
----i
----Loop Bound: (1, 1024)
------j
------Loop Bound: (0, 1024)
--------array access ey.addr ((i * 1024) + j)
--------array access hz.addr ((i * 1024) + j)
--------array access hz.addr (((i - 1) * 1024) + j)
--------array access ey.addr ((i * 1024) + j)
----i
----Loop Bound: (0, 1024)
------j
------Loop Bound: (1, 1024)
--------array access ex.addr ((i * 1024) + j)
--------array access hz.addr ((i * 1024) + j)
--------array access hz.addr (((i * 1024) + j) - 1)
--------array access ex.addr ((i * 1024) + j)
----i
----Loop Bound: (0, 1023)
------j
------Loop Bound: (0, 1023)
--------array access hz.addr ((i * 1024) + j)
--------array access ex.addr (((i * 1024) + j) + 1)
--------array access ex.addr ((i * 1024) + j)
--------array access ey.addr (((i + 1) * 1024) + j)
--------array access ey.addr ((i * 1024) + j)
--------array access hz.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
fdtd_2d */ 
