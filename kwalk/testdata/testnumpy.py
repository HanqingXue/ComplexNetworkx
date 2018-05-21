import numpy as np

a = []
a.append([1, 1, 1])
a.append([2, 1, 2])
a.append([2, 2, 2])
print a
a = np.array([, , ])
print np.apply_along_axis(sum, 1, a)
print a.sum(axis=1)