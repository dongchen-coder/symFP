
 /* Start to analysis array index
Array index info
E.addr ((i * 256) + j)
A.addr ((i * 256) + k)
B.addr ((k * 256) + j)
E.addr ((i * 256) + j)
E.addr ((i * 256) + j)
F.addr ((i * 256) + j)
C.addr ((i * 256) + k)
D.addr ((k * 256) + j)
F.addr ((i * 256) + j)
F.addr ((i * 256) + j)
G.addr ((i * 256) + j)
E.addr ((i * 256) + k)
F.addr ((k * 256) + j)
G.addr ((i * 256) + j)
G.addr ((i * 256) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %ni
i32 %nj
i32 %nk
i32 %nl
i32 %nm
double* %E
double* %A
double* %B
double* %F
double* %C
double* %D
double* %G

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 256)
----j
----Loop Bound: (0, 256)
------array access E.addr ((i * 256) + j)
------k
------Loop Bound: (0, 256)
--------array access A.addr ((i * 256) + k)
--------array access B.addr ((k * 256) + j)
--------array access E.addr ((i * 256) + j)
--------array access E.addr ((i * 256) + j)
--i
--Loop Bound: (0, 256)
----j
----Loop Bound: (0, 256)
------array access F.addr ((i * 256) + j)
------k
------Loop Bound: (0, 256)
--------array access C.addr ((i * 256) + k)
--------array access D.addr ((k * 256) + j)
--------array access F.addr ((i * 256) + j)
--------array access F.addr ((i * 256) + j)
--i
--Loop Bound: (0, 256)
----j
----Loop Bound: (0, 256)
------array access G.addr ((i * 256) + j)
------k
------Loop Bound: (0, 256)
--------array access E.addr ((i * 256) + k)
--------array access F.addr ((k * 256) + j)
--------array access G.addr ((i * 256) + j)
--------array access G.addr ((i * 256) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
mm2 */ 
