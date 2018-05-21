import numpy as np  
import time
import gnumpy as gpu  

def bench_np():
	n = 4000  
	for i in range(10):  
	    a = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)    
	    b = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)  
	    a = a.dot(b)  

def bench_gnp():
	n = 4000  
	for i in range(10):  
	    a = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)    
	    b = np.random.uniform(low=0., high=1., size=(n, n)).astype(np.float32)  
	    ga = gpu.garray(a)  
	    gb = gpu.garray(b)  
	  
	    ga = ga.dot(gb)

def main():
	t1 = time.time()
	bench_np()
	t2 = time.time()


	t3 = time.time()
	bench_gnp()
	t4 = time.time()

	print 'Numpy:{}'.format(t2 - t1)
	print 'Gumpy:{}'.format(t4 - t3)

if __name__ == '__main__':
  	main()  