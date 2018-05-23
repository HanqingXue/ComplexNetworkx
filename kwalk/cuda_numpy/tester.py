import numpy as np
import gnumpy as gpu 
from timeit import timeit 

def test_numpy():
    n = 10000  
    for i in range(10):  
        a = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)    
        b = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)  
        a = a.dot(b)  


def test_gnumpy():
    n = 10000
    for i in range(10):  
        a = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)    
        b = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)  
        ga = gpu.garray(a)  
        gb = gpu.garray(b)  
      
        ga = ga.dot(gb)

cuda_time = timeit('test_gnumpy()', 'from __main__ import test_gnumpy', number=1)
cpu_time = timeit('test_numpy()', 'from __main__ import test_numpy', number=1)

print cpu_time
print cuda_time

