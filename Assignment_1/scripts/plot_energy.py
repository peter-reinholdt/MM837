#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

x = np.loadtxt(filename)

E = x[:,1]

plt.figure()
plt.plot(E - E[0])
plt.xlabel("Timestep (a.u.)")
plt.ylabel(r"$\Delta E$")
plt.savefig("energy_error.png")

plt.figure()

plt.plot(E[:100] - E[0])
plt.xlabel("Timestep (a.u.)")
plt.ylabel(r"$\Delta E$")
plt.savefig("energy_error_start.png")
