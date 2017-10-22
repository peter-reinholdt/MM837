#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

timesteps = np.logspace(-4,0,200)
errors    = np.zeros(200)

for i in range(200):
    print(i)
    x = np.loadtxt("{}.dat".format(i))
    errors[i] = np.max(np.abs(x[:,1] - x[0,1]))

plt.plot(timesteps, errors)
plt.xscale("log")
plt.yscale("log")
plt.show()
