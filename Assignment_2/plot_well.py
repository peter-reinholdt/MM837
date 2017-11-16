#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np

x = np.loadtxt(sys.argv[1], delimiter=",")
plt.plot(x[:,0], x[:,1], label=r"$\psi$")
plt.plot(x[:,0], x[:,2], label=r"$V$")
plt.plot(x[:,0], x[:,3], label=r"$E$")

plt.xlabel(r"$\tilde{x}$")
plt.ylabel("Property")
plt.legend()
plt.show()
