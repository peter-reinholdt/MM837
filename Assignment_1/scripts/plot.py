#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

x = np.loadtxt("output.dat")

plt.plot(x[:,1])
plt.show()
