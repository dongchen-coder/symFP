
 /* Start to analysis array index
Array index info
data.addr ((i * 1024) + j)
mean.addr j
data.addr ((i * 1024) + j)
mean.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
stddev.addr j
mean.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
stddev.addr j
data.addr ((i * 1024) + j)
data.addr ((i * 1024) + j)
symmat.addr ((j1 * 1024) + j1)
symmat.addr ((j1 * 1024) + j2)
data.addr ((i * 1024) + j1)
data.addr ((i * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j1 * 1024) + j2)
symmat.addr ((j2 * 1024) + j1)
symmat.addr 1048575
mean.addr j
data.addr ((i * 1024) + j)
mean.addr j
mean.addr j
mean.addr j
mean.addr j
stddev.addr j

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %m
i32 %n
double* %data
double* %mean
double* %stddev
double* %symmat

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--j
--Loop Bound: (0, 1024)
----array access mean.addr j
----i
----Loop Bound: (0, 1024)
------array access data.addr ((i * 1024) + j)
------array access mean.addr j
------array access mean.addr j
----array access mean.addr j
----array access mean.addr j
--j
--Loop Bound: (0, 1024)
----array access stddev.addr j
----i
----Loop Bound: (0, 1024)
------array access data.addr ((i * 1024) + j)
------array access mean.addr j
------array access data.addr ((i * 1024) + j)
------array access mean.addr j
------array access stddev.addr j
------array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
----array access stddev.addr j
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access mean.addr j
------array access data.addr ((i * 1024) + j)
------array access data.addr ((i * 1024) + j)
------array access stddev.addr j
------array access data.addr ((i * 1024) + j)
------array access data.addr ((i * 1024) + j)
--j1
--Loop Bound: (0, 1023)
----array access symmat.addr ((j1 * 1024) + j1)
----j2
----Loop Bound: ((j1 + 1), 1024)
------array access symmat.addr ((j1 * 1024) + j2)
------i
------Loop Bound: (0, 1024)
--------array access data.addr ((i * 1024) + j1)
--------array access data.addr ((i * 1024) + j2)
--------array access symmat.addr ((j1 * 1024) + j2)
--------array access symmat.addr ((j1 * 1024) + j2)
------array access symmat.addr ((j1 * 1024) + j2)
------array access symmat.addr ((j2 * 1024) + j1)
--array access symmat.addr 1048575

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
correlation */ 
