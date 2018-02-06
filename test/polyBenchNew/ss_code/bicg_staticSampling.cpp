
 /* Start to analysis array index
Array index info
q.addr i
s.addr i
s.addr j
r.addr i
A.addr ((i * 1024) + j)
s.addr j
q.addr i
A.addr ((i * 1024) + j)
p.addr j
q.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %nx
i32 %ny
double* %A
double* %r
double* %s
double* %p
double* %q

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----array access s.addr i
--i
--Loop Bound: (0, 1024)
----array access q.addr i
----j
----Loop Bound: (0, 1024)
------array access s.addr j
------array access r.addr i
------array access A.addr ((i * 1024) + j)
------array access s.addr j
------array access q.addr i
------array access A.addr ((i * 1024) + j)
------array access p.addr j
------array access q.addr i

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
bicg_cpu */ 
