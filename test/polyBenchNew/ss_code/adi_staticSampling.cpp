
 /* Start to analysis array index
Array index info
v.addr (0 + i)
p.addr ((i * 1024) + 0)
v.addr (0 + i)
q.addr ((i * 1024) + 0)
p.addr (((i * 1024) + j) - 1)
p.addr ((i * 1024) + j)
u.addr (((j * 1024) + i) - 1)
u.addr ((j * 1024) + i)
u.addr (((j * 1024) + i) + 1)
q.addr (((i * 1024) + j) - 1)
p.addr (((i * 1024) + j) - 1)
q.addr ((i * 1024) + j)
v.addr (1047552 + i)
p.addr ((i * 1024) + j)
v.addr (((j + 1) * 1024) + i)
q.addr ((i * 1024) + j)
v.addr ((j * 1024) + i)
u.addr ((i * 1024) + 0)
p.addr ((i * 1024) + 0)
u.addr ((i * 1024) + 0)
q.addr ((i * 1024) + 0)
p.addr (((i * 1024) + j) - 1)
p.addr ((i * 1024) + j)
v.addr (((i - 1) * 1024) + j)
v.addr ((i * 1024) + j)
v.addr (((i + 1) * 1024) + j)
q.addr (((i * 1024) + j) - 1)
p.addr (((i * 1024) + j) - 1)
q.addr ((i * 1024) + j)
u.addr (((i * 1024) + 1024) - 1)
p.addr ((i * 1024) + j)
u.addr (((i * 1024) + j) + 1)
q.addr ((i * 1024) + j)
u.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %p
double* %q
double* %v
double* %u

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (1, 10)
----i
----Loop Bound: (1, 1023)
------array access v.addr (0 + i)
------array access p.addr ((i * 1024) + 0)
------array access v.addr (0 + i)
------array access q.addr ((i * 1024) + 0)
------j
------Loop Bound: (1, 1023)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access p.addr ((i * 1024) + j)
--------array access u.addr (((j * 1024) + i) - 1)
--------array access u.addr ((j * 1024) + i)
--------array access u.addr (((j * 1024) + i) + 1)
--------array access q.addr (((i * 1024) + j) - 1)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access q.addr ((i * 1024) + j)
------array access v.addr (1047552 + i)
------j
------Loop Bound: (1022, 1)
--------array access p.addr ((i * 1024) + j)
--------array access v.addr (((j + 1) * 1024) + i)
--------array access q.addr ((i * 1024) + j)
--------array access v.addr ((j * 1024) + i)
----i
----Loop Bound: (1, 1023)
------array access u.addr ((i * 1024) + 0)
------array access p.addr ((i * 1024) + 0)
------array access u.addr ((i * 1024) + 0)
------array access q.addr ((i * 1024) + 0)
------j
------Loop Bound: (1, 1023)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access p.addr ((i * 1024) + j)
--------array access v.addr (((i - 1) * 1024) + j)
--------array access v.addr ((i * 1024) + j)
--------array access v.addr (((i + 1) * 1024) + j)
--------array access q.addr (((i * 1024) + j) - 1)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access q.addr ((i * 1024) + j)
------array access u.addr (((i * 1024) + 1024) - 1)
------j
------Loop Bound: (1022, 1)
--------array access p.addr ((i * 1024) + j)
--------array access u.addr (((i * 1024) + j) + 1)
--------array access q.addr ((i * 1024) + j)
--------array access u.addr ((i * 1024) + j)

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
 /* Start to analyze function:  
adi */ 
