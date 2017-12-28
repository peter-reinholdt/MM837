#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

filenames = sys.argv[1:]

data = np.zeros(np.loadtxt(filenames[0]).shape)


plt.subplot(211)
for fn in filenames:
    x = np.loadtxt(fn)
    data += x
    plt.plot(x[:,0], x[:,1], alpha=0.1)
data = data / len(filenames)
plt.plot(data[:,0], data[:,1])
plt.subplot(212)

for fn in filenames:
    x = np.loadtxt(fn)
    plt.plot(x[:,0], x[:,2], alpha=0.1)

plt.plot(data[:,0], data[:,2])
plt.show()
