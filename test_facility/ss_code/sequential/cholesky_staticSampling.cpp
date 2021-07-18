
 /* Start to analysis array index
Array index info: Total number of references: 13
A.addr ((i * 1024) + k)
A.addr ((j * 1024) + k)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + i)
A.addr ((j * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + i)
A.addr ((i * 1024) + i)
A.addr ((i * 1024) + i)
BC Array cost info: Total number of arrays: 1
A.addr 24
BC Average cost: 2.400000e+01

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
Array A_addr 
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, i)
----Loop inc: (j + 1)
----Loop predicate: <
------k
------Loop Bound: (0, j)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A_addr ((i * 1024) + k)
--------array access A_addr ((j * 1024) + k)
--------array access A_addr ((i * 1024) + j)
--------array access A_addr ((i * 1024) + j)
------array access A_addr ((j * 1024) + j)
------array access A_addr ((i * 1024) + j)
------array access A_addr ((i * 1024) + j)
----k
----Loop Bound: (0, i)
----Loop inc: (k + 1)
----Loop predicate: <
------array access A_addr ((i * 1024) + k)
------array access A_addr ((i * 1024) + k)
------array access A_addr ((i * 1024) + i)
------array access A_addr ((i * 1024) + i)
----array access A_addr ((i * 1024) + i)
----array access A_addr ((i * 1024) + i)

Finish analysis loops */ 
/* Start IV Dependence Analysis i:
  - j(2)
  - k(2)
j:
  - k(3)
Finish to analyze IV Dependence *//* Start to analysis the access graphA_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr
A_addr -> A_addr

*/
 /* Start to analysis the number of samples
calculating:
init counter: 0 0 
Dump stride: 1 1 
init counter: 0 0 0 
Dump stride: 1 1 1 
