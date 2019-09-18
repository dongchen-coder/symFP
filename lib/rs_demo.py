import numpy as np

# An efficient Python3 program  
# to randomly select k items 
# from a stream of items 
import random 
# A utility function  
# to print an array 
def printArray(stream,n): 
    for i in range(n): 
        print(stream[i],end=" "); 
    print(); 
  
# A function to randomly select 
# k items from stream[0..n-1]. 
def selectKItems(stream, n, k): 
        i=0;  
        # index for elements 
        # in stream[] 
          
        # reservoir[] is the output  
        # array. Initialize it with 
        # first k elements from stream[] 
        reservoir = [0]*k; 
        for i in range(k): 
            reservoir[i] = stream[i]; 
          
        # Iterate from the (k+1)th 
        # element to nth element 
        while(i < n): 
            # Pick a random index 
            # from 0 to i. 
            j = random.randrange(i+1); 
              
            # If the randomly picked 
            # index is smaller than k, 
            # then replace the element 
            # present at the index 
            # with new element from stream 
            if(j < k): 
                reservoir[j] = stream[i]; 
            i+=1; 
          
        print("Following are k randomly selected items"); 
        printArray(reservoir, k); 


# n:    number of resverior
def kernel(n, dim_size):
    # init the counter for each resverior
    counter = [0] * n
    reservoirs = [0] * n
    reservoir_cnt = n
    for i in xrange(0, dim_size + 1):
        for j in xrange(0, dim_size + 1):
            b[i* (dim_size + 2) +j] =  a[i* (dim_size + 2)+j] 
                                    + a[i* (dim_size + 2)+j + 1] 
                                    + a[i* (dim_size + 2)+j - 1] 
                                    + a[(i-1)* (dim_size + 2) +j] 
                                    + a[(i+1)* (dim_size + 2) +j];
            # init the resveriors
            if reservoir_cnt == n:
                reservoirs[n - reservoir_cnt] = i * (dim_size + 2) + j
                reservoir_cnt -= 1
                reservoirs[n - reservoir_cnt] = i * (dim_size + 2) + j + 1
                reservoir_cnt -= 1
                reservoirs[n - reservoir_cnt] = i * (dim_size + 2) + j - 1
                reservoir_cnt -= 1
                reservoirs[n - reservoir_cnt] = (i-1) * (dim_size + 2) + j
                reservoir_cnt -= 1
                reservoirs[n - reservoir_cnt] = (i+1) * (dim_size + 2) + j
                reservoir_cnt -= 1
                reservoirs[n - reservoir_cnt] = i * (dim_size + 2) + j + (dim_size + 2)* (dim_size + 2)
                reservoir_cnt -= 1
            else:
                for r in xrange(0, n):
                    if i * (dim_size + 2) + j

