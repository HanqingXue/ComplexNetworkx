import numpy as np  
n = 4000  
for i in range(10):  
    a = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)    b = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)  
    a = a.dot(b)  