#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]
x = np.loadtxt(filename)
plt.plot(x[:,0], x[:,1])
plt.savefig(filename + ".png")
