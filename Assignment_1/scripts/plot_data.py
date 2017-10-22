#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]

stop = 20000

x = np.loadtxt(filename)
plt.figure(figsize=(10,5))
plt.subplot(421)
plt.plot(x[:stop,0], "k.")
plt.xlabel("Timestep (x100)")
plt.ylabel(r"$<p^2>$")



plt.subplot(422)
plt.plot(x[:stop,1])
plt.xlabel("Timestep (x100)")
plt.ylabel(r"$E$")

plt.subplot(423)
plt.plot(x[:stop,2])
plt.xlabel("Timestep (x100)")
plt.ylabel(r"$T$")


plt.subplot(424)
plt.plot(x[:stop,3])
plt.xlabel("Timestep (x100)")
plt.ylabel(r"$U$")


plt.subplot(425)
plt.plot(x[:stop,4])
plt.xlabel("Timestep (x100)")
plt.ylabel(r"$P$")


plt.subplot(426)
plt.plot(x[:stop,5])
plt.xlabel("Timestep (x100)")
plt.ylabel(r"$Q$")



plt.subplot(427)
plt.plot(x[:10000,0])
plt.xlabel("Timestep (x100)")
plt.ylabel(r"$<p^2>$")


plt.subplot(428)
plt.plot(x[:1000,0])
plt.xlabel("Timestep (x100)")
plt.ylabel(r"$<p^2>$")




plt.figure()
plt.hist(x[10000:,0], bins=100)
plt.xlabel(r"$<p^2>$")
plt.ylabel("Count")
plt.show()

print(np.average(x[10000:,0]), np.var(x[10000:,0]))
