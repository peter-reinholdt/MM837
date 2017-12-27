#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]

data = np.loadtxt(filename)

plt.subplot(211)
plt.plot(data[:,0], data[:,1])
plt.subplot(212)
plt.plot(data[:,0], data[:,2])
plt.show()

